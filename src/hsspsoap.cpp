//  Copyright Maarten L. Hekkelman, Radboud University 2008.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

//#define NDEBUG
#undef DEBUG
#define DEBUG 0

#include "MRS.h"

#include <sys/wait.h>
#include <unistd.h>
#include <pwd.h>

#include <boost/algorithm/string.hpp>
#include <boost/filesystem/operations.hpp>
#include <boost/filesystem/convenience.hpp>
#include <boost/filesystem/fstream.hpp>
#include <boost/format.hpp>
#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/filter/zlib.hpp>
#include <boost/iostreams/device/back_inserter.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/copy.hpp>
#include <boost/program_options.hpp>
#include <boost/lexical_cast.hpp>

#include "zeep/config.hpp"

#include "CDatabank.h"
#include "CDatabankTable.h"
#include "CBlast.h"
#include "CQuery.h"

#include "zeep/server.hpp"

#include "blast.h"
#include "maxhom-hssp.h"

#define HSSPSOAP_PID_FILE	"/var/run/hsspsoap.pid"
#define HSSPSOAP_LOG_FILE	"/var/log/hsspsoap.log"

using namespace std;
namespace ba = boost::algorithm;
namespace io = boost::iostreams;
namespace po = boost::program_options;

// globals

fs::path gMaxHom;

void GetDSSPForPDBFile(
	const string&		inPDBFile,
	string&				outDSSP)
{
	HUuid uuid;
	
	fs::path rundir("/tmp/hsspsoap/");
	rundir /= boost::lexical_cast<string>(uuid);
	fs::create_directories(rundir);

	// write out the dssp file
	fs::ofstream pdb(rundir / "in.pdb");
	pdb << inPDBFile;
	pdb.close();
	
	fs::path dsspcmbi = "/home/maarten/projects/dssp/dsspcmbi";
	if (not fs::exists(dsspcmbi))
		THROW(("dsspcmbi executable '%s' not found", dsspcmbi.string().c_str()));
	
	int pid = fork();
	if (pid == -1)
		THROW(("fork failed: %s", strerror(errno)));
	
	if (pid == 0)	// the child process (will be maxhom)
	{
		fs::current_path(rundir);
		
		setpgid(0, 0);
		
		int fd = open("dsspcmbi.log", O_CREAT | O_RDWR, 0666);
		dup2(fd, STDOUT_FILENO);
		dup2(fd, STDERR_FILENO);
		
		char* argv[] = {
			const_cast<char*>(dsspcmbi.string().c_str()),
			const_cast<char*>("in.pdb"),
			const_cast<char*>("out.dssp"),
			nil
		};
		
		(void)execve(dsspcmbi.string().c_str(), argv, environ);
	}

	int status;
	waitpid(pid, &status, 0);
	
	if (status != 0)
		THROW(("dsspcmbi exited with status %d", status));
	
	// OK, got it! Read in the result and exit
	fs::ifstream result(rundir / "out.dssp");
	io::filtering_ostream out(io::back_inserter(outDSSP));
	io::copy(result, out);
}

void GetDSSPForSequence(
	const string&		inSequence,
	string&				outDSSP)
{
	stringstream s;

	char first[] =	"==== Secondary Structure Definition by the program DSSP, updated CMBI version by ElmK / April 1,2000 ==== DATE=28-MAY-2010     .\n"
					"REFERENCE\n"
					"HEADER\n"
					"COMPND\n"
					"SOURCE\n"
					"AUTHOR\n"
					"    x  1  0  0  0\n"
					"  #\n";
	
	char *p = first + sizeof(first) - 19;
	if (*p != 'x')
		THROW(("Fout! p = '%s'", p));

	for (int n = inSequence.length(); n > 0; n /= 10)
		*p-- = '0' + (n % 10);
	
	s << first;
	
	// And now fill in the rest
	
	int n = 1;
	for (string::const_iterator aa = inSequence.begin(); aa != inSequence.end(); ++aa, ++n)
	{
		char line[17] = {};
		
		snprintf(line, sizeof(line), "%5.d%5.d A %c", n, n, toupper(*aa));
		s << line << endl;
	}

	outDSSP = s.str();
}

void ParseSequenceFromDSSP(
	const string&		inDSSP,
	string&				outSequence)
{
	istringstream s(inDSSP);
	
	string line;
	getline(s, line);
	if (not ba::starts_with(line, "==== Secondary Structure Definition"))
		THROW(("Not a valid DSSP file"));
	
	int n = 0, state = 0;
	
	outSequence.clear();
	outSequence.reserve(n);
	
	while (state < 100)
	{
		getline(s, line);
		if (s.eof() or line.empty())
			THROW(("Invalid (truncated?) DSSP file"));
		
		switch (state)
		{
			case 0:
				if (ba::starts_with(line, "REFERENCE"))
				{
					state = 1;
					break;
				}

			case 1:
				if (ba::starts_with(line, "HEADER"))
				{
					state = 2;
					break;
				}

			case 2:
				if (ba::starts_with(line, "COMPND"))
				{
					state = 3;
					break;
				}

			case 3:
				if (ba::starts_with(line, "SOURCE"))
				{
					state = 4;
					break;
				}
			
			case 4:
				if (ba::starts_with(line, "AUTHOR"))
				{
					state = 5;
					break;
				}

			case 5:
				if (ba::starts_with(line, " "))
				{
					line.erase(5, string::npos);
					while (not line.empty() and isspace(line[0]))
						line.erase(line.begin());
					n = boost::lexical_cast<int>(line);
					state = 6;
				}
				break;
			
			case 6:
				if (ba::starts_with(line, "  #"))
					state = 100;
				break;
		}
	}
	
	outSequence.reserve(n);
	while (n-- > 0 and not s.eof())
	{
		getline(s, line);
		
		if (line.length() < 14)
			THROW(("Invalid DSSP file"));
		char aa = line[13];
		if (islower(aa))
			aa = 'C';
		outSequence.push_back(aa);
	}
}

// the data types used in our communication with the outside world
// are wrapped in a namespace.

class hssp_server : public zeep::server
{
  public:
			hssp_server();

	void	GetDSSPForPDBID(
				const string&	pdbid,
				string&			dssp);

	void	GetDSSPForPDBFile(
				const string&	pdbfile,
				string&			dssp);

	void	GetHSSPForPDBID(
				const string&	pdbid,
				string&			hssp);

	void	GetHSSPForPDBFile(
				const string&	pdbfile,
				string&			hssp);

	void	GetHSSPForSequence(
				const string&	sequence,
				string&			hssp);

	CDatabankTable				mDBTable;
};

hssp_server::hssp_server()
	: zeep::server("http://www.cmbi.ru.nl/hsspsoap", "hsspsoap")
{
//	const char* kGetDSSPForPDBIDParameterNames[] = {
//		"pdbid", "dssp"
//	};
//	
//	register_action("GetDSSPForPDBID", this, &hssp_server::GetDSSPForPDBID, kGetDSSPForPDBIDParameterNames);

	const char* kGetDSSPForPDBFileParameterNames[] = {
		"pdbfile", "dssp"
	};
	
	register_action("GetDSSPForPDBFile", this, &hssp_server::GetDSSPForPDBFile, kGetDSSPForPDBFileParameterNames);

//	const char* kGetHSSPForPDBIDParameterNames[] = {
//		"pdbid", "hssp"
//	};
//	
//	register_action("GetHSSPForPDBID", this, &hssp_server::GetHSSPForPDBID, kGetHSSPForPDBIDParameterNames);

	const char* kGetHSSPForPDBFileParameterNames[] = {
		"pdbfile", "hssp"
	};
	
	register_action("GetHSSPForPDBFile", this, &hssp_server::GetHSSPForPDBFile, kGetHSSPForPDBFileParameterNames);

	const char* kGetHSSPForSequenceParameterNames[] = {
		"sequence", "hssp"
	};
	
	register_action("GetHSSPForSequence", this, &hssp_server::GetHSSPForSequence, kGetHSSPForSequenceParameterNames);
}

//void hssp_server::GetDSSPForPDBID(
//	const string&				pdbid,
//	string&						dssp)
//{
//	string sequence;
//	
//	GetSequenceForPDBID(pdbid, sequence);
//	
//	vector<uint32> hits;
//	BlastSequence(sequence, hits);
//	
//	dssp = boost::lexical_cast<string>(hits.size());
//}

void hssp_server::GetDSSPForPDBFile(
	const string&				pdbfile,
	string&						dssp)
{
	::GetDSSPForPDBFile(pdbfile, dssp);
}

//void hssp_server::GetHSSPForPDBID(
//	const string&				pdbid,
//	string&						hssp)
//{
//}

void hssp_server::GetHSSPForPDBFile(
	const string&				pdbfile,
	string&						hssp)
{
	CDatabankPtr db = mDBTable.Load("uniprot");

	string dssp;
	::GetDSSPForPDBFile(pdbfile, dssp);
	
	string sequence;
	::ParseSequenceFromDSSP(dssp, sequence);
	
	vector<uint32> hits;
	::BlastSequence(db, sequence, hits);
	
	ostringstream s;
	maxhom::GetHSSPForHitsAndDSSP(db, gMaxHom.string(), "UNKN", hits, dssp, 1500, s);
	hssp = s.str();
}

void hssp_server::GetHSSPForSequence(
	const string&				sequence,
	string&						hssp)
{
	CDatabankPtr db = mDBTable.Load("uniprot");

	string dssp;
	::GetDSSPForSequence(sequence, dssp);
	
	vector<uint32> hits;
	BlastSequence(db, sequence, hits);
	
	ostringstream s;
	maxhom::GetHSSPForHitsAndDSSP(db, gMaxHom.string(), "UNKN", hits, dssp, 1500, s);
	hssp = s.str();
}

// --------------------------------------------------------------------
//
//	Daemonize
// 

void Daemonize(
	const string&		inUser)
{
	int pid = fork();
	
	if (pid == -1)
	{
		cerr << "Fork failed" << endl;
		exit(1);
	}
	
	if (pid != 0)
		_exit(0);

	if (setsid() < 0)
	{
		cerr << "Failed to create process group: " << strerror(errno) << endl;
		exit(1);
	}

	// it is dubious if this is needed:
	signal(SIGHUP, SIG_IGN);

	// fork again, to avoid being able to attach to a terminal device
	pid = fork();

	if (pid == -1)
		cerr << "Fork failed" << endl;

	if (pid != 0)
		_exit(0);

	// write our pid to the pid file
	ofstream pidFile(HSSPSOAP_PID_FILE);
	pidFile << getpid() << endl;
	pidFile.close();

	if (chdir("/") != 0)
	{
		cerr << "Cannot chdir to /: " << strerror(errno) << endl;
		exit(1);
	}

	if (inUser.length() > 0)
	{
		struct passwd* pw = getpwnam(inUser.c_str());
		if (pw == NULL or setuid(pw->pw_uid) < 0)
		{
			cerr << "Failed to set uid to " << inUser << ": " << strerror(errno) << endl;
			exit(1);
		}
	}

	// close stdin
	close(STDIN_FILENO);
	open("/dev/null", O_RDONLY);
}

// --------------------------------------------------------------------
// 
//	OpenLogFile
//	

void OpenLogFile()
{
	// open the log file
	int fd = open(HSSPSOAP_LOG_FILE, O_CREAT|O_APPEND|O_RDWR, 0644);
	if (fd < 0)
	{
		cerr << "Opening log file " HSSPSOAP_LOG_FILE " failed" << endl;
		exit(1);
	}

	// redirect stdout and stderr to the log file
	dup2(fd, STDOUT_FILENO);
	dup2(fd, STDERR_FILENO);
}

// --------------------------------------------------------------------
// 
//	main
//	

int main(int argc, char* argv[])
{
	po::options_description desc("Options");
	desc.add_options()
		("help,h",								"Display help message")
		("address,a",	po::value<string>(),	"address to bind to")
		("port,p",		po::value<uint16>(),	"port to bind to")
		("location,l",	po::value<string>(),	"location advertised in wsdl")
		("user,u",		po::value<string>(),	"user to run as")
		("maxhom",		po::value<string>(),	"Path to the maxhom application")
		("threads,a",	po::value<int>(),		"Number of threads to use (default is nr of CPU's)")
		("no-daemon,D",							"do not fork a daemon")
		;
	
	string
		location = "http://mrs.cmbi.ru.nl/hsspsoap/wsdl",
		address = "0.0.0.0",
		user = "nobody";

	uint16
		port = 10334;

	bool
		daemon = true;

	po::variables_map vm;
	po::store(po::parse_command_line(argc, argv, desc), vm);
	po::notify(vm);
	
	if (vm.count("help"))
	{
		cout << desc << endl;
		exit(1);
	}
	
	if (vm.count("address"))
		address = vm["address"].as<string>();

	if (vm.count("location"))
		location = vm["location"].as<string>();

	if (vm.count("port"))
		port = vm["port"].as<uint16>();

	if (vm.count("user"))
		user = vm["user"].as<string>();

	string maxhom = "maxhom";
	if (vm.count("maxhom"))
		maxhom = vm["maxhom"].as<string>();

	if (vm.count("threads"))
		BLAST_THREADS = vm["threads"].as<int>();
	
	gMaxHom = maxhom;
	if (not fs::exists(gMaxHom))
	{
		cerr << "No maxhom found" << endl;
		exit(1);
	}

	if (vm.count("no-daemon"))
		daemon = false;

	if (daemon)
	{
		OpenLogFile();
		Daemonize(user);
	}
	
#ifndef _MSC_VER
    sigset_t new_mask, old_mask;
    sigfillset(&new_mask);
    pthread_sigmask(SIG_BLOCK, &new_mask, &old_mask);
#endif

	hssp_server server;
	server.bind(address, port);
	
	if (not location.empty())
		server.set_location(location);
	
    boost::thread t(boost::bind(&hssp_server::run, &server, 1));

#ifndef _MSC_VER
    pthread_sigmask(SIG_SETMASK, &old_mask, 0);

	// Wait for signal indicating time to shut down.
	sigset_t wait_mask;
	sigemptyset(&wait_mask);
	sigaddset(&wait_mask, SIGINT);
	sigaddset(&wait_mask, SIGQUIT);
	sigaddset(&wait_mask, SIGTERM);
	pthread_sigmask(SIG_BLOCK, &wait_mask, 0);
	int sig = 0;
	sigwait(&wait_mask, &sig);
	
	server.stop();
#endif

	t.join();

	return 0;
}
