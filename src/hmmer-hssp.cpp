//  Copyright Maarten L. Hekkelman, Radboud University 2011.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include "MRS.h"

#if P_UNIX
#include <wait.h>
#elif P_WIN
#include <Windows.h>
#endif

#include <cmath>

#include <boost/algorithm/string.hpp>
#include <boost/filesystem/convenience.hpp>
#include <boost/filesystem/fstream.hpp>
#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/copy.hpp>
#include <boost/iostreams/filter/bzip2.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/foreach.hpp>
#define foreach BOOST_FOREACH
#include <boost/date_time/gregorian/gregorian.hpp>
#include <boost/date_time/date_clock_device.hpp>
#include <boost/regex.hpp>
#include <boost/ptr_container/ptr_vector.hpp>

// MRS includes
#include "CDatabank.h"
#include "CUtils.h"
#include "CConfig.h"

// our includes
#include "matrix.h"
#include "dssp.h"
#include "structure.h"
#include "utils.h"
#include "hmmer-hssp.h"

using namespace std;
namespace ba = boost::algorithm;
namespace io = boost::iostreams;

namespace hmmer
{

// global, 5 minutes
uint32 gMaxRunTime = 300, gNrOfThreads = boost::thread::hardware_concurrency();

// precalculated threshold table for identity values between 10 and 80
const double kHomologyThreshold[] = {
	0.845468, 0.80398,  0.767997, 0.736414, 0.708413, 0.683373, 0.660811, 0.640351, 0.621688, 0.604579,
	0.58882,  0.574246, 0.560718, 0.548117, 0.536344, 0.525314, 0.514951, 0.505194, 0.495984, 0.487275,
	0.479023, 0.471189, 0.463741, 0.456647, 0.449882, 0.44342,  0.43724,  0.431323, 0.425651, 0.420207,
	0.414976, 0.409947, 0.405105, 0.40044,  0.395941, 0.391599, 0.387406, 0.383352, 0.379431, 0.375636,
	0.37196,  0.368396, 0.364941, 0.361587, 0.358331, 0.355168, 0.352093, 0.349103, 0.346194, 0.343362,
	0.340604, 0.337917, 0.335298, 0.332744, 0.330252, 0.327821, 0.325448, 0.323129, 0.320865, 0.318652,
	0.316488, 0.314372, 0.312302, 0.310277, 0.308294, 0.306353, 0.304452, 0.302589, 0.300764, 0.298975,
	0.297221,
};

// --------------------------------------------------------------------
// utility routine
	
inline bool is_gap(char aa)
{
	return aa == '-' or aa == '~' or aa == '.' or aa == '_';
}

void SetMaxRunTime(uint32 inSeconds)
{
	gMaxRunTime = inSeconds;
}

void SetNrOfThreads(uint32 inThreads)
{
	gNrOfThreads = inThreads;
}

// --------------------------------------------------------------------
// basic named sequence type and a multiple sequence alignment container
	
struct seq
{
	string		m_id;
	string		m_seq;
	uint32		m_identical, m_length;
	uint32		m_begin, m_end;
	
				seq() {}

				seq(const string& id)
					: m_id(id)
					, m_identical(0)
					, m_length(0)
					, m_begin(0)
					, m_end(0)
				{
					m_seq.reserve(5000);
				}
};

typedef boost::ptr_vector<seq> mseq;

// --------------------------------------------------------------------
// ReadStockholm is a function that reads a multiple sequence alignment from
// a Stockholm formatted file. Restriction is that this Stockholm file has
// a #=GF field at the second line containing the ID of the query used in
// jackhmmer.	

void ReadStockholm(istream& is, mseq& msa)
{
	if (VERBOSE)
		cerr << "Reading stockholm file...";

	string line, qseq;
	getline(is, line);
	if (line != "# STOCKHOLM 1.0")
		throw mas_exception("Not a stockholm file");

	getline(is, line);
	if (not ba::starts_with(line, "#=GF ID "))
		throw mas_exception("Not a valid stockholm file, missing #=GF ID line");
	
	string id = line.substr(8);

	boost::regex re("(.+?)-i(?:\\d+)$");
	boost::smatch sm;
	if (boost::regex_match(id, sm, re))
		id = sm.str(1);

	msa.push_back(new seq(id));
	uint32 ix = 0;
	
	for (;;)
	{
		getline(is, line);
		
		if (line.empty())
		{
			if (is.eof())
				break;
			continue;
		}
		
		if (line == "//")
			break;
		
		if (ba::starts_with(line, "#=GS "))
		{
			
			string id = line.substr(5);
			string::size_type s = id.find("DE ");
			if (s != string::npos)
				id = id.substr(0, s);
			
			ba::trim(id);
			if (msa.size() > 1 or msa.front().m_id != id)
				msa.push_back(new seq(id));
			continue;
		}
		
		if (line[0] != '#')
		{
			string::size_type s = line.find(' ');
			if (s == string::npos)
				throw mas_exception("Invalid stockholm file");
			
			string id = line.substr(0, s);
			
			while (s < line.length() and line[s] == ' ')
				++s;
			
			string seq = line.substr(s);
			
			if (id == msa[0].m_id)
			{
				ix = 0;
				msa[0].m_seq += seq;
				qseq = seq;
			}
			else
			{
				++ix;
				if (ix >= msa.size())
					msa.push_back(new hmmer::seq(id));

				assert(ix < msa.size());
				assert(id == msa[ix].m_id);
				
				msa[ix].m_seq += seq;
				
				for (string::const_iterator qi = qseq.begin(), si = seq.begin(); qi != qseq.end(); ++qi, ++si)
				{
					if (not is_gap(*qi) and *qi == *si)
						++msa[ix].m_identical;
					
					if (not is_gap(*qi) or not is_gap(*si))
						++msa[ix].m_length;
				}
			}
		}
	}
	
	if (msa.size() < 2)
		THROW(("Insufficient sequences in Stockholm MSA"));

	if (VERBOSE)
		cerr << " done" << endl << "Checking for threshold...";
	
	// for our query
	string& q = msa.front().m_seq;
	msa.front().m_begin = 0;
	msa.front().m_end = q.length();

	// Remove all hits that are not above the threshold here
	mseq::iterator mi = msa.begin() + 1;
	while (mi != msa.end())
	{
		float score = float(mi->m_identical) / float(mi->m_length);
		uint32 ix = max(10U, min(mi->m_length, 80U)) - 10;
		if (score < kHomologyThreshold[ix])
		{
			if (VERBOSE > 1)
				cerr << "dropping " << mi->m_id << " because identity " << score << " is below threshold " << kHomologyThreshold[ix] << endl;
			
			mi = msa.erase(mi);
		}
		else
		{
			string& s = mi->m_seq;
			
			// update begin/end here
			for (uint32 i = 0; i < s.length(); ++i)
			{
				if (is_gap(s[i]))
					continue;

				if (mi->m_begin == 0 and i != 0)
					mi->m_begin = i;
				
				mi->m_end = i + 1;
				assert(mi->m_end <= s.length());
				assert(mi->m_end <= q.length());
			}

			assert(mi->m_end <= s.length());
			assert(mi->m_end <= q.length());

			++mi;
		}
	}

	if (VERBOSE)
		cerr << "done" << endl;
}

void CheckAlignmentForChain(
	mseq&			inMSA,
	const MChain*	inChain)
{
	string sa, sc;

	foreach (char r, inMSA.front().m_seq)
	{
		if (not is_gap(r))
			sa += r;
	}

	inChain->GetSequence(sc);

	if (sa != sc)
	{
		if (sa.length() < sc.length())
			THROW(("Query used for Stockholm file is too short for the chain"));

		string::size_type offset = sa.find(sc);
		if (offset == string::npos)
			THROW(("Invalid Stockholm file for chain"));

		if (offset > 0)
		{
			foreach (seq& s, inMSA)
			{
				s.m_seq.erase(0, offset);
				if (s.m_begin > offset)
					s.m_begin -= offset;
				else
					s.m_begin = 0;
				
				if (s.m_end > offset)
					s.m_end -= offset;
				else
					s.m_end = 0;
			}
		}

		if (sa.length() > sc.length() + offset)
		{
			uint32 n = sa.length() - (sc.length() + offset);
			foreach (seq& s, inMSA)
			{
				uint32 o = s.m_seq.length() - n;
				
				s.m_seq.erase(o, n);
				
				if (s.m_begin > o)
					s.m_begin = o;
				if (s.m_end > o)
					s.m_end = o;
			}
		}
	}
}

#if P_UNIX
// --------------------------------------------------------------------
// Run the Jackhmmer application

fs::path RunJackHmmer(const string& seq, uint32 iterations, const fs::path& fastadir,
	const fs::path& jackhmmer, const string& db)
{
	if (seq.empty())
		THROW(("Empty sequence in RunJackHmmer"));
	
	HUuid uuid;
	
	fs::path rundir("/tmp/hssp-2/");
	rundir /= boost::lexical_cast<string>(uuid);
	fs::create_directories(rundir);
	
	if (VERBOSE)
		cerr << "Running jackhmmer (" << uuid << ")...";
		
	// write fasta file
	fs::ofstream input(rundir / "input.fa");
	if (not input.is_open())
		throw mas_exception("Failed to create jackhmmer input file");
		
	input << '>' << "input" << endl;
	for (uint32 o = 0; o < seq.length(); o += 72)
	{
		uint32 k = seq.length() - o;
		if (k > 72)
			k = 72;
		input << seq.substr(o, k) << endl;
	}
	input.close();
	
	// start a jackhmmer
	int pid = fork();
	
	if (pid == -1)
		THROW(("fork failed: %s", strerror(errno)));
	
	if (pid == 0)	// the child process (will be jackhmmer)
	{
		fs::current_path(rundir);
		
		setpgid(0, 0);
		
		int fd = open("jackhmmer.log", O_CREAT | O_RDWR | O_APPEND, 0666);
		dup2(fd, STDOUT_FILENO);
		dup2(fd, STDERR_FILENO);
		
		arg_vector argv(jackhmmer.string());
		
		argv.push("-N", iterations);
		argv.push("--noali");
		argv.push("--cpu", gNrOfThreads);
//		argv.push("-o", "/dev/null");
		argv.push("-A", "output.sto");
		argv.push("input.fa");
		argv.push((fastadir / (db + ".fa")).string());
		
		if (VERBOSE)
			cerr << argv << endl;
		
		(void)execve(jackhmmer.string().c_str(), argv, environ);
		cerr << "Failed to run " << jackhmmer << endl << " err: " << strerror(errno) << endl;
		exit(-1);
	}

	// wait for jackhmmer to finish or time out
	double startTime = system_time();
	int status;

	for (;;)
	{
		int err = waitpid(pid, &status, WNOHANG);
		if (err == -1 or err == pid)
			break;
		
		if (system_time() > startTime + gMaxRunTime)
		{
			err = kill(pid, SIGKILL);
			if (err == 0)
				err = waitpid(pid, &status, 0);
			
			THROW(("Timeout waiting for jackhmmer result"));
		}
		
		sleep(1);
	}
	
	if (status != 0)
	{
		if (fs::exists(rundir / "jackhmmer.log"))
		{
			fs::ifstream log(rundir / "jackhmmer.log");
			
			if (log.is_open())
			{
				// only print the last 10 lines
				deque<string> lines;
			
				for (;;)
				{
					string line;
					getline(log, line);
					
					if (line.empty() and log.eof())
						break;
					
					lines.push_back(line);
					if (lines.size() > 10)
						lines.pop_front();
				}
				
				foreach (string& line, lines)
					cerr << line;
			}
		}
		
		THROW(("jackhmmer exited with status %d", status));
	}

	if (not fs::exists(rundir / "output.sto"))
		THROW(("Output Stockholm file is missing"));

	return rundir;
}

#elif P_WIN

fs::path RunJackHmmer(const string& seq, uint32 iterations, const fs::path& fastadir,
	const fs::path& jackhmmer, const string& db)
{
	// Jackhmmer as downloaded from http://hmmer.janelia.org/software is a cygwin application
	// this means we can use 

	if (seq.empty())
		THROW(("Empty sequence in RunJackHmmer"));
	
	HUuid uuid;
	
	fs::path rundir(gScratchDir / "hssp-2");
	rundir /= boost::lexical_cast<string>(uuid);
	fs::create_directories(rundir);
	
	if (VERBOSE)
		cerr << "Running jackhmmer (" << uuid << ")...";
		
	// write fasta file
	
	fs::ofstream input(rundir / "input.fa");
	if (not input.is_open())
		throw mas_exception("Failed to create jackhmmer input file");
		
	input << '>' << "input" << endl;
	for (uint32 o = 0; o < seq.length(); o += 72)
	{
		uint32 k = seq.length() - o;
		if (k > 72)
			k = 72;
		input << seq.substr(o, k) << endl;
	}
	input.close();
	
	static int sSerial = 1;

	// fork/exec a jackhmmer to do the work
	if (not fs::exists(jackhmmer))
		THROW(("The jackhmmer executable '%s' does not seem to exist", jackhmmer.string().c_str()));

	double startTime = system_time();
	
	// ready to roll
	SECURITY_ATTRIBUTES sa = { sizeof(SECURITY_ATTRIBUTES) };
	sa.bInheritHandle = true;

	enum { i_read, i_write, i_write2 };
	HANDLE ifd[3], ofd[3], efd[3];

	::CreatePipe(&ifd[i_read], &ifd[i_write], &sa, 0);
	::DuplicateHandle(::GetCurrentProcess(), ifd[i_write],
		::GetCurrentProcess(), &ifd[i_write2], 0, false,
		DUPLICATE_SAME_ACCESS);
	::CloseHandle(ifd[i_write]);

	::CreatePipe(&ofd[i_read], &ofd[i_write], &sa, 0);
	::DuplicateHandle(::GetCurrentProcess(), ofd[i_read],
		::GetCurrentProcess(), &ofd[i_write2], 0, false,
		DUPLICATE_SAME_ACCESS);
	::CloseHandle(ofd[i_read]);

	::CreatePipe(&efd[i_read], &efd[i_write], &sa, 0);
	::DuplicateHandle(::GetCurrentProcess(), efd[i_read],
		::GetCurrentProcess(), &efd[i_write2], 0, false,
		DUPLICATE_SAME_ACCESS);
	::CloseHandle(efd[i_read]);

	STARTUPINFOA si = { sizeof(STARTUPINFOA) };
	si.dwFlags = STARTF_USESHOWWINDOW | STARTF_USESTDHANDLES;
	si.hStdInput = ifd[i_read];
	si.hStdOutput = ofd[i_write];
	si.hStdError = efd[i_write];

	string cwd = rundir.string();

	stringstream scmd;
	scmd << jackhmmer << ' '
		<< "-N " << iterations << ' '
		<< "--noali" << ' '
		<< "--cpu " << gNrOfThreads << ' '
		<< "-A " << (rundir / "output.sto") << ' '
		<< (rundir / "input.fa") << ' '
		<< (fastadir / (db + ".fa"));
	string cmd = scmd.str();

	PROCESS_INFORMATION pi;
	::CreateProcessA(nil, const_cast<char*>(cmd.c_str()), nil, nil, true,
		CREATE_NEW_PROCESS_GROUP, nil, const_cast<char*>(cwd.c_str()), &si, &pi);

	::CloseHandle(ifd[i_read]);
	::CloseHandle(ofd[i_write]);
	::CloseHandle(efd[i_write]);

	HANDLE proc = pi.hProcess;
	HANDLE thread = pi.hThread;
	DWORD pid = pi.dwProcessId;
	DWORD tid = pi.dwThreadId;

	DWORD rr, avail;

	// OK, so now the executable is started and the pipes are set up
	// write the sequences and read from the pipes until done.
	bool errDone = false, outDone = false, killed = false;
	string error, out;
	
	while (not errDone and not outDone)
	{
		::Sleep(100);

		char buffer[1024];

		while (not outDone)
		{
			if (not ::PeekNamedPipe(ofd[i_write2], nil, 0, nil, &avail, nil))
			{
				unsigned int err = ::GetLastError();
				if (err == ERROR_HANDLE_EOF or err == ERROR_BROKEN_PIPE)
					outDone = true;
			}
			else if (avail > 0 and ::ReadFile(ofd[i_write2], buffer, sizeof(buffer), &rr, nil))
				out.append(buffer, buffer + rr);
			else
				break;
		}

		while (not errDone)
		{
			if (not ::PeekNamedPipe(efd[i_write2], nil, 0, nil, &avail, nil))
			{
				unsigned int err = ::GetLastError();
				if (err == ERROR_HANDLE_EOF or err == ERROR_BROKEN_PIPE)
					errDone = true;
			}
			else if (avail > 0 and ::ReadFile(efd[i_write2], buffer, sizeof(buffer), &rr, nil))
				error.append(buffer, buffer + rr);
			else
				break;
		}

		if (not errDone and not outDone and not killed and startTime + gMaxRunTime < system_time())
		{
			::TerminateProcess(proc, 1);

			// is this enough?
			::CloseHandle(ofd[i_write2]);
			::CloseHandle(efd[i_write2]);
			::CloseHandle(proc);
			::CloseHandle(thread);

			THROW(("jackhmmer was killed since its runtime exceeded the limit of %d seconds", gMaxRunTime));
		}
	}

	::CloseHandle(ofd[i_write2]);
	::CloseHandle(efd[i_write2]);
	::CloseHandle(proc);
	::CloseHandle(thread);

	if (not error.empty())
		cerr << error << endl;

	if (not fs::exists(rundir / "output.sto"))
		THROW(("Output Stockholm file is missing"));

	return rundir;
}

#endif

void RunJackHmmer(const string& seq, uint32 iterations, const fs::path& fastadir, const fs::path& jackhmmer,
	const string& db, fs::path dst)
{
	fs::path rundir = RunJackHmmer(seq, iterations, fastadir, jackhmmer, db);

	// copy the result

	fs::ifstream in(rundir / "output.sto");
	fs::ofstream outfile(dst, ios_base::binary);

	io::filtering_stream<io::output> out;
	if (dst.extension() == ".bz2")
		out.push(io::bzip2_compressor());
	else if (dst.extension() == ".gz")
		out.push(io::gzip_compressor());
	out.push(outfile);

	io::copy(in, out);

	if (not VERBOSE)
		fs::remove_all(rundir);
	else
		cerr << " done" << endl;
}

void RunJackHmmer(const string& seq, uint32 iterations, const fs::path& fastadir, const fs::path& jackhmmer,
	const string& db, mseq& msa)
{
	fs::path rundir = RunJackHmmer(seq, iterations, fastadir, jackhmmer, db);

	fs::ifstream is(rundir / "output.sto");
	ReadStockholm(is, msa);
	is.close();

	// read in the result
	if (not fs::exists(rundir / "output.sto"))
		THROW(("Output Stockholm file is missing"));
	
	if (not VERBOSE)
		fs::remove_all(rundir);
	else
		cerr << " done" << endl;
}

// --------------------------------------------------------------------
// Hit is a class to store hit information and all of its statistics.
	
struct insertion
{
	uint32			ipos, jpos;
	string			seq;
};
	
struct Hit
{
					Hit(mseq& msa, char chain, uint32 qix, uint32 six);

	const mseq&		msa;
	char			chain;
	uint32			nr, ix;
	string			id, acc, desc, pdb;
	uint32			ifir, ilas, jfir, jlas, lali, ngap, lgap, lseq2;
	float			ide, wsim;
	uint32			identical, similar;

	bool			operator<(const Hit& rhs) const 	{ return ide > rhs.ide or (ide == rhs.ide and lali > rhs.lali); }
	
	vector<insertion>
					insertions;
};

typedef shared_ptr<Hit> hit_ptr;

// Create a Hit object based on a jackhmmer alignment pair
// first is the original query sequence, with gaps introduced.
// second is the hit sequence.
// Since this is jackhmmer output, we can safely assume the
// alignment does not contain gaps at the start or end of the query.
Hit::Hit(mseq& msa, char chain, uint32 qix, uint32 six)
	: msa(msa)
	, chain(chain)
	, ix(six)
{
	string& q = msa[qix].m_seq;
	string& s = msa[six].m_seq;

	assert(q.length() == s.length());

	uint32 b = msa[six].m_begin;
	uint32 e = msa[six].m_end;

	assert(b < q.length());
	assert(e <= q.length());

	assert(not q.empty() and not s.empty());
	assert(not is_gap(q[0]) and not is_gap(q[q.length() - 1]));
	
	// parse out the position
	static const boost::regex re("([-a-zA-Z0-9_]+)/(\\d+)-(\\d+)");
	boost::smatch sm;

	if (not boost::regex_match(msa[six].m_id, sm, re))
		throw mas_exception("Alignment ID should contain position");
	
	id = sm.str(1);
	
	ifir = b + 1;
	ilas = b;

	// jfir/jlas can be taken over from jackhmmer output
	jfir = boost::lexical_cast<uint32>(sm.str(2));
	jlas = boost::lexical_cast<uint32>(sm.str(3));

	lgap = ngap = identical = similar = 0;
	lali = s.length();
	
	string::iterator qb = q.begin() + b, qe = q.begin() + e,
					 sb = s.begin() + b, se = s.begin() + e;

	fill(s.begin(), sb, ' ');
	fill(se, s.end(), ' ');

	bool sgap = false, qgap = false;
	const substitution_matrix m("BLOSUM62");
	uint32 ipos = ifir, jpos = jfir;
	insertion ins;
	
	for (string::iterator si = sb, qi = qb; si != se; ++si, ++qi)
	{
		if (is_gap(*si) and is_gap(*qi))
		{
			// a common gap
			--lali;
		}
		else if (is_gap(*si))
		{
			if (not (sgap or qgap))
				++ngap;
			sgap = true;
			++ilas;
			++lgap;
			++jpos;
		}
		else if (is_gap(*qi))
		{
			if (not qgap)
			{
				assert(si != sb);
				string::iterator gsi = si - 1;
				while (gsi != sb and is_gap(*gsi))
					--gsi;
				
				*gsi = tolower(*gsi);
				ins.ipos = ipos;
				ins.jpos = jpos - 1;
				ins.seq.assign(gsi, gsi + 1);
				ins.seq += *si;
			}
			else
				ins.seq += *si;
			
			if (not (sgap or qgap))
				++ngap;

			qgap = true;
			++lgap;
			++ipos;
		}
		else
		{
			if (qgap)
			{
				*si = tolower(*si);
				ins.seq += *si;
				insertions.push_back(ins);
			}
			
			sgap = false;
			qgap = false;

			if (*qi == *si)
			{
				++identical;
				++similar;
			}
			else if (m(*qi, *si) > 0)
				++similar;
			
			++ilas;
			++ipos;
			++jpos;
		}
	}

	ide = float(identical) / float(lali);
	wsim = float(similar) / float(lali);
}

struct compare_hit
{
	bool operator()(hit_ptr a, hit_ptr b) const { return *a < *b; }
};

// --------------------------------------------------------------------
// ResidueHInfo is a class to store information about a residue in the
// original query sequence, along with statistics.

struct ResidueHInfo
{
					ResidueHInfo(uint32 seqNr);
					ResidueHInfo(char a, vector<hit_ptr>& hits, uint32 pos, char chain, uint32 seqNr, uint32 pdbNr,
						const string& dssp, float consweight);
	
	char			letter;
	char			chain;
	string			dssp;
	uint32			seqNr, pdbNr;
	uint32			pos;
	uint32			nocc, ndel, nins;
	float			entropy, consweight;
	uint32			dist[20];
	static int8		kIX[256];
};

typedef shared_ptr<ResidueHInfo>	res_ptr;

// --------------------------------------------------------------------

//const string ResidueHInfo::kIX("VLIMFWYGAPSTCHRKQEND");

int8 ResidueHInfo::kIX[256] = {
	//   1   2   3   4   5   6   7   8   9  10  11  12  13  14  15
	-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, //  0
	-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, //  1
	-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, //  2 
	-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, //  3 
	-1,  8, -1, 12, 19, 17,  4,  7, 13,  2, -1, 15,  1,  3, 18, -1, //  4 
	 9, 16, 14, 10, 11, -1,  0,  5, -1,  6, -1, -1, -1, -1, -1, -1, //  5 
	-1,  8, -1, 12, 19, 17,  4,  7, 13,  2, -1, 15,  1,  3, 18, -1, //  4 
	 9, 16, 14, 10, 11, -1,  0,  5, -1,  6, -1, -1, -1, -1, -1, -1, //  5 
	-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, //  8 
	-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, //  9 
	-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, // 10 
	-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, // 11 
	-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, // 12 
	-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, // 13 
	-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, // 14 
	-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1  // 15 
};


// first constructor is for a 'chain-break'
ResidueHInfo::ResidueHInfo(uint32 seqNr)
	: letter(0)
	, seqNr(seqNr)
{
}

ResidueHInfo::ResidueHInfo(char a, vector<hit_ptr>& hits, uint32 pos, char chain, uint32 seqNr, uint32 pdbNr,
		const string& dssp, float consweight)
	: letter(a)
	, chain(chain)
	, dssp(dssp)
	, seqNr(seqNr)
	, pdbNr(pdbNr)
	, pos(pos)
	, nocc(1)
	, ndel(0)
	, nins(0)
	, consweight(consweight)
{
	fill(dist, dist + 20, 0);
	
	int8 ix = kIX[uint8(a)];
	assert(ix != -1);
	if (ix != -1)
		dist[ix] = 1;
	
	foreach (hit_ptr hit, hits)
	{
		ix = kIX[uint8(hit->msa[hit->ix].m_seq[pos])];
		if (ix != -1)
		{
			++nocc;
			dist[ix] += 1;
		}
	}
	
	entropy = 0;
	for (uint32 a = 0; a < 20; ++a)
	{
		double freq = double(dist[a]) / nocc;
		
		dist[a] = uint32((100.0 * freq) + 0.5);
		
		if (freq > 0)
			entropy -= static_cast<float>(freq * log(freq));
	}
	
	// calculate ndel and nins
	const mseq& msa = hits.front()->msa;
	const string& q = msa[0].m_seq;
	
	bool gap = pos + 1 < q.length() and is_gap(q[pos + 1]);
	
	foreach (hit_ptr hit, hits)
	{
		const string& t = msa[hit->ix].m_seq;
		
		if (is_gap(t[pos]))
			++ndel;
		
		if (gap and t[pos] >= 'a' and t[pos] <= 'y')
			++nins;
	}
}

// --------------------------------------------------------------------
// Write collected information as a HSSP file to the output stream

void CreateHSSPOutput(
	const string&		inProteinID,
	const string&		inProteinDescription,
	const string&		inDatabankVersion,
	uint32				inSeqLength,
	uint32				inNChain,
	uint32				inKChain,
	const string&		inUsedChains,
	vector<hit_ptr>&	hits,
	vector<res_ptr>&	res,
	ostream&			os)
{
	using namespace boost::gregorian;
	date today = day_clock::local_day();
	
	// print the header
	os << "HSSP       HOMOLOGY DERIVED SECONDARY STRUCTURE OF PROTEINS , VERSION 2.0d2 2011" << endl
	   << "PDBID      " << inProteinID << endl
	   << "DATE       file generated on " << to_iso_extended_string(today) << endl
	   << "SEQBASE    " << inDatabankVersion << endl
	   << "THRESHOLD  according to: t(L)=(290.15 * L ** -0.562) + 5" << endl
	   << "CONTACT    This version: Maarten L. Hekkelman <m.hekkelman@cmbi.ru.nl>" << endl
	   << inProteinDescription
	   << boost::format("SEQLENGTH  %4.4d") % inSeqLength << endl
	   << boost::format("NCHAIN     %4.4d chain(s) in %s data set") % inNChain % inProteinID << endl;
	
	if (inKChain != inNChain)
	{
		os << boost::format("KCHAIN     %4.4d chain(s) used here ; chains(s) : ") % inKChain << inUsedChains << endl;
	}
	
	os << boost::format("NALIGN     %4.4d") % hits.size() << endl
	   << endl
	   << "## PROTEINS : identifier and alignment statistics" << endl
	   << "  NR.    ID         STRID   %IDE %WSIM IFIR ILAS JFIR JLAS LALI NGAP LGAP LSEQ2 ACCNUM     PROTEIN" << endl;
	   
	// print the first list
	uint32 nr = 1;
	boost::format fmt1("%5.5d : %12.12s%4.4s    %4.2f  %4.2f %4.4d %4.4d %4.4d %4.4d %4.4d %4.4d %4.4d %4.4d  %10.10s %s");
	foreach (hit_ptr h, hits)
	{
		string id = h->id;
		if (id.length() > 12)
			id.erase(12, string::npos);
		else if (id.length() < 12)
			id.append(12 - id.length(), ' ');
		
		string acc = h->acc;
		if (acc.length() > 10)
			acc.erase(10, string::npos);
		else if (acc.length() < 10)
			acc.append(10 - acc.length(), ' ');

		os << fmt1 % nr
				   % id % h->pdb
				   % h->ide % h->wsim % h->ifir % h->ilas % h->jfir % h->jlas % h->lali
				   % h->ngap % h->lgap % h->lseq2
				   % acc % h->desc
		   << endl;
		
		++nr;
	}

	// print the alignments
	for (uint32 i = 0; i < hits.size(); i += 70)
	{
		uint32 n = i + 70;
		if (n > hits.size())
			n = hits.size();
		
		uint32 k[7] = {
			((i +  0) / 10) % 10 + 1,
			((i + 10) / 10) % 10 + 1,
			((i + 20) / 10) % 10 + 1,
			((i + 30) / 10) % 10 + 1,
			((i + 40) / 10) % 10 + 1,
			((i + 50) / 10) % 10 + 1,
			((i + 60) / 10) % 10 + 1
		};
		
		os << boost::format("## ALIGNMENTS %4.4d - %4.4d") % (i + 1) % n << endl
		   << boost::format(" SeqNo  PDBNo AA STRUCTURE BP1 BP2  ACC NOCC  VAR  ....:....%1.1d....:....%1.1d....:....%1.1d....:....%1.1d....:....%1.1d....:....%1.1d....:....%1.1d")
		   					% k[0] % k[1] % k[2] % k[3] % k[4] % k[5] % k[6] << endl;

		res_ptr last;
		foreach (res_ptr ri, res)
		{
			if (ri->letter == 0)
				os << boost::format(" %5.5d        !  !           0   0    0    0    0") % ri->seqNr << endl;
			else
			{
				string aln;
				
				for (uint32 j = i; j < n; ++j)
				{
					if (hits[j]->chain == ri->chain)
					{
						uint32 p = ri->pos;
						uint32 i = hits[j]->ix;
						aln += hits[j]->msa[i].m_seq[p];
					}
					else
						aln += ' ';
				}
				
				uint32 ivar = uint32(100 * (1 - ri->consweight));

				os << ' ' << boost::format("%5.5d%s%4.4d %4.4d  ") % ri->seqNr % ri->dssp % ri->nocc % ivar << aln << endl;
			}
		}
	}
	
	// ## SEQUENCE PROFILE AND ENTROPY
	os << "## SEQUENCE PROFILE AND ENTROPY" << endl
	   << " SeqNo PDBNo   V   L   I   M   F   W   Y   G   A   P   S   T   C   H   R   K   Q   E   N   D  NOCC NDEL NINS ENTROPY RELENT WEIGHT" << endl;
	
	res_ptr last;
	foreach (res_ptr r, res)
	{
		if (r->letter == 0)
		{
			os << boost::format("%5.5d          0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0     0    0    0   0.000      0")
				% r->seqNr << endl;
		}
		else
		{
			os << boost::format(" %4.4d %4.4d %c") % r->seqNr % r->pdbNr % r->chain;

			for (uint32 i = 0; i < 20; ++i)
				os << boost::format("%4.4d") % r->dist[i];

			uint32 relent = uint32(100 * r->entropy / log(20.0));
			os << "  " << boost::format("%4.4d %4.4d %4.4d   %5.3f   %4.4d  %4.2f") % r->nocc % r->ndel % r->nins % r->entropy % relent % r->consweight << endl;
		}
	}
	
	// insertion list
	
	os << "## INSERTION LIST" << endl
	   << " AliNo  IPOS  JPOS   Len Sequence" << endl;

	foreach (hit_ptr h, hits)
	{
		foreach (insertion& ins, h->insertions)
		{
			string s = ins.seq;
			
			if (s.length() <= 100)
				os << boost::format("  %4.4d  %4.4d  %4.4d  %4.4d ") % h->nr % ins.ipos % ins.jpos % (ins.seq.length() - 2) << s << endl;
			else
			{
				os << boost::format("  %4.4d  %4.4d  %4.4d  %4.4d ") % h->nr % ins.ipos % ins.jpos % (ins.seq.length() - 2) << s.substr(0, 100) << endl;
				s.erase(0, 100);
				
				while (not s.empty())
				{
					uint32 n = s.length();
					if (n > 100)
						n = 100;
					
					os << "     +                   " << s.substr(0, n) << endl;
					s.erase(0, n);
				}
			}
		}			
	}
	
	os << "//" << endl;
}

// --------------------------------------------------------------------
// Calculate the variability of a residue, based on dayhoff similarity
// and weights

// Dayhoff matrix as used by maxhom
const float kDayhoffData[] =
{
     1.5f,                                                                                                                  // V
     0.8f, 1.5f,                                                                                                            // L
     1.1f, 0.8f, 1.5f,                                                                                                      // I
     0.6f, 1.3f, 0.6f, 1.5f,                                                                                                // M
     0.2f, 1.2f, 0.7f, 0.5f, 1.5f,                                                                                          // F
    -0.8f, 0.5f,-0.5f,-0.3f, 1.3f, 1.5f,                                                                                    // W
    -0.1f, 0.3f, 0.1f,-0.1f, 1.4f, 1.1f, 1.5f,                                                                              // Y
     0.2f,-0.5f,-0.3f,-0.3f,-0.6f,-1.0f,-0.7f, 1.5f,                                                                        // G
     0.2f,-0.1f, 0.0f, 0.0f,-0.5f,-0.8f,-0.3f, 0.7f, 1.5f,                                                                  // A
     0.1f,-0.3f,-0.2f,-0.2f,-0.7f,-0.8f,-0.8f, 0.3f, 0.5f, 1.5f,                                                            // P
    -0.1f,-0.4f,-0.1f,-0.3f,-0.3f, 0.3f,-0.4f, 0.6f, 0.4f, 0.4f, 1.5f,                                                      // S
     0.2f,-0.1f, 0.2f, 0.0f,-0.3f,-0.6f,-0.3f, 0.4f, 0.4f, 0.3f, 0.3f, 1.5f,                                                // T
     0.2f,-0.8f, 0.2f,-0.6f,-0.1f,-1.2f, 1.0f, 0.2f, 0.3f, 0.1f, 0.7f, 0.2f, 1.5f,                                          // C
    -0.3f,-0.2f,-0.3f,-0.3f,-0.1f,-0.1f, 0.3f,-0.2f,-0.1f, 0.2f,-0.2f,-0.1f,-0.1f, 1.5f,                                    // H
    -0.3f,-0.4f,-0.3f, 0.2f,-0.5f, 1.4f,-0.6f,-0.3f,-0.3f, 0.3f, 0.1f,-0.1f,-0.3f, 0.5f, 1.5f,                              // R
    -0.2f,-0.3f,-0.2f, 0.2f,-0.7f, 0.1f,-0.6f,-0.1f, 0.0f, 0.1f, 0.2f, 0.2f,-0.6f, 0.1f, 0.8f, 1.5f,                        // K
    -0.2f,-0.1f,-0.3f, 0.0f,-0.8f,-0.5f,-0.6f, 0.2f, 0.2f, 0.3f,-0.1f,-0.1f,-0.6f, 0.7f, 0.4f, 0.4f, 1.5f,                  // Q
    -0.2f,-0.3f,-0.2f,-0.2f,-0.7f,-1.1f,-0.5f, 0.5f, 0.3f, 0.1f, 0.2f, 0.2f,-0.6f, 0.4f, 0.0f, 0.3f, 0.7f, 1.5f,            // E
    -0.3f,-0.4f,-0.3f,-0.3f,-0.5f,-0.3f,-0.1f, 0.4f, 0.2f, 0.0f, 0.3f, 0.2f,-0.3f, 0.5f, 0.1f, 0.4f, 0.4f, 0.5f, 1.5f,      // N
    -0.2f,-0.5f,-0.2f,-0.4f,-1.0f,-1.1f,-0.5f, 0.7f, 0.3f, 0.1f, 0.2f, 0.2f,-0.5f, 0.4f, 0.0f, 0.3f, 0.7f, 1.0f, 0.7f, 1.5f // D
};

// --------------------------------------------------------------------
// Convert a multiple sequence alignment as created by jackhmmer to 
// a set of information as used by HSSP.

void ChainToHits(CDatabankPtr inDatabank, mseq& msa, const MChain& chain,
	vector<hit_ptr>& hits, vector<res_ptr>& res)
{
	if (VERBOSE)
		cerr << "Creating hits...";
	
	vector<hit_ptr> nhits;

	for (uint32 i = 1; i < msa.size(); ++i)
	{
		hit_ptr h(new Hit(msa, chain.GetChainID(), 0, i));

		uint32 docNr = inDatabank->GetDocumentNr(h->id);
		
		h->desc = inDatabank->GetMetaData(docNr, "title");
		try
		{
			if (ba::starts_with(h->id, "UniRef100_"))
				h->acc = h->id = h->id.substr(10);
			else
				h->acc = inDatabank->GetMetaData(docNr, "acc");
		}
		catch (...) {}
		h->lseq2 = inDatabank->GetSequence(docNr, 0).length();
		
		// update number now that we know how far we are
		h->ifir += res.size();
		h->ilas += res.size();
		
		nhits.push_back(h);
	}
	
	if (VERBOSE)
		cerr << " done" << endl
			 << "Continuing with " << nhits.size() << " hits" << endl
			 << "Calculating conservation weights...";

	const string& s = msa.front().m_seq;
	vector<float> sumvar(s.length()), sumdist(s.length()), simval(s.length());
	static const symmetric_matrix<float> D(kDayhoffData, 20);

	for (uint32 i = 0; i + 1 < msa.size(); ++i)
	{
		for (uint32 j = i + 1; j < msa.size(); ++j)
		{
			const string& si = msa[i].m_seq;
			const string& sj = msa[j].m_seq;

			uint32 b = msa[i].m_begin;
			if (b < msa[j].m_begin)
				b = msa[j].m_begin;
			
			uint32 e = msa[i].m_end;
			if (e > msa[j].m_end)
				e = msa[j].m_end;

			uint32 len = 0, agr = 0;
			for (uint32 k = b; k < e; ++k)
			{
				if (not is_gap(si[k]) and not is_gap(sj[k]))
				{
					++len;
					if (si[k] == sj[k])
						++agr;

					int8 ri = ResidueHInfo::kIX[uint8(si[k])];
					int8 rj = ResidueHInfo::kIX[uint8(sj[k])];
					
					if (ri != -1 and rj != -1)
						simval[k] = D(ri, rj);
					else
						simval[k] = numeric_limits<float>::min();
				}
			}

			if (len > 0)
			{
				float distance = 1 - (float(agr) / float(len));
				for (uint32 k = b; k < e; ++k)
				{
					if (simval[k] != numeric_limits<float>::min())
					{
						sumvar[k] += distance * simval[k];
						sumdist[k] += distance * 1.5f;
					}
				}
			}
		}
	}

	if (VERBOSE)
		cerr << " done" << endl
			 << "Calculating residue info...";

	const vector<MResidue*>& residues = chain.GetResidues();
	vector<MResidue*>::const_iterator ri = residues.begin();

	for (uint32 i = 0; i < s.length(); ++i)
	{
		if (is_gap(s[i]))
			continue;

		assert(ri != residues.end());
		
		if (ri != residues.begin() and (*ri)->GetNumber() > (*(ri - 1))->GetNumber() + 1)
			res.push_back(res_ptr(new ResidueHInfo(res.size() + 1)));
		
		string dssp = ResidueToDSSPLine(**ri).substr(5, 34);

		float weight = 1.0f;
		if (sumdist[i] > 0)
			weight = sumvar[i] / sumdist[i];

		res.push_back(res_ptr(new ResidueHInfo(s[i], nhits, i,
			chain.GetChainID(), res.size() + 1, (*ri)->GetNumber(), dssp, weight)));

		++ri;
	}
	
	if (VERBOSE)
		cerr << " done" << endl;
	
	assert(ri == residues.end());
	hits.insert(hits.end(), nhits.begin(), nhits.end());
}

// Find the minimal set of overlapping sequences
// Only search fully contained subsequences, no idea what to do with
// sequences that overlap and each have a tail. What residue number to use in that case? What chain ID?
void ClusterSequences(vector<string>& s, vector<uint32>& ix)
{
	for (;;)
	{
		bool found = false;
		for (uint32 i = 0; not found and i < s.size() - 1; ++i)
		{
			for (uint32 j = i + 1; not found and j < s.size(); ++j)
			{
				string& a = s[i];
				string& b = s[j];

				if (a.empty() or b.empty())
					continue;

				if (ba::contains(a, b)) // j fully contained in i
				{
					s[j].clear();
					ix[j] = i;
					found = true;
				}
				else if (ba::contains(b, a)) // i fully contained in j
				{
					s[i].clear();
					ix[i] = j;
					found = true;
				}
			}
		}
		
		if (not found)
			break;
	}
}

void CreateHSSP(
	CDatabankPtr				inDatabank,
	const string&				inProtein,
	const fs::path&				inFastaDir,
	const fs::path&				inJackHmmer,
	uint32						inIterations,
	ostream&					outHSSP)
{
	vector<hit_ptr> hits;
	vector<res_ptr> res;
	mseq alignment;

	RunJackHmmer(inProtein, inIterations, inFastaDir, inJackHmmer, inDatabank->GetID(), alignment);

	MChain chain('A');
	vector<MResidue*>& residues = chain.GetResidues();
	MResidue* last = nil;
	uint32 nr = 1;
	foreach (char r, inProtein)
	{
		residues.push_back(new MResidue(nr, r, last));
		++nr;
		last = residues.back();
	}
	
	ChainToHits(inDatabank, alignment, chain, hits, res);

	sort(hits.begin(), hits.end(), compare_hit());
	if (hits.size() > 9999)
		hits.erase(hits.begin() + 9999, hits.end());
	
	nr = 1;
	foreach (hit_ptr h, hits)
		h->nr = nr++;

	CreateHSSPOutput("UNKN", "", inDatabank->GetVersion(), inProtein.length(),
		1, 1, "A", hits, res, outHSSP);
}

void CreateHSSP(
	CDatabankPtr				inDatabank,
	MProtein&					inProtein,
	const fs::path&				inFastaDir,
	const fs::path&				inJackHmmer,
	uint32						inIterations,
	uint32						inMinSeqLength,
	ostream&					outHSSP)
{
	uint32 seqlength = 0;

	vector<hit_ptr> hits;
	vector<res_ptr> res;

	// construct a set of unique sequences, containing only the largest ones in case of overlap
	vector<string> seqset;
	vector<uint32> ix;
	vector<const MChain*> chains;
	uint32 kchain = 0;
	
	foreach (const MChain* chain, inProtein.GetChains())
	{
		string seq;
		chain->GetSequence(seq);
		
		if (seq.length() < inMinSeqLength)
			continue;
		
		chains.push_back(chain);
		seqset.push_back(seq);
		ix.push_back(ix.size());
	}
	
	if (seqset.empty())
		THROW(("Not enough sequences in DSSP file of length %d", inMinSeqLength));

	if (seqset.size() > 1)
		ClusterSequences(seqset, ix);
	
	// only take the unique sequences
	ix.erase(unique(ix.begin(), ix.end()), ix.end());

	// Maybe we should change this code to run jackhmmer only once 
	vector<mseq> alignments(seqset.size());
	foreach (uint32 i, ix)
		RunJackHmmer(seqset[i], inIterations, inFastaDir, inJackHmmer, inDatabank->GetID(), alignments[i]);
	
	foreach (uint32 i, ix)
	{
		const MChain* chain = chains[i];
		
		string& seq = seqset[i];
		assert(not seq.empty());
		seqlength += seq.length();

		if (not res.empty())
			res.push_back(res_ptr(new ResidueHInfo(res.size() + 1)));

		ChainToHits(inDatabank, alignments[i], *chain, hits, res);
		++kchain;
	}

	sort(hits.begin(), hits.end(), compare_hit());
	if (hits.size() > 9999)
		hits.erase(hits.begin() + 9999, hits.end());
	
	uint32 nr = 1;
	foreach (hit_ptr h, hits)
		h->nr = nr++;
	
	string usedChains;
	foreach (uint32 i, ix)
	{
		if (not usedChains.empty())
			usedChains += ',';
		usedChains += chains[i]->GetChainID();
	}
	
	stringstream desc;
	desc	
	   << "HEADER     " + inProtein.GetHeader().substr(10, 40) << endl
	   << "COMPND     " + inProtein.GetCompound().substr(10) << endl
	   << "SOURCE     " + inProtein.GetSource().substr(10) << endl
	   << "AUTHOR     " + inProtein.GetAuthor().substr(10) << endl;

	CreateHSSPOutput(inProtein.GetID(), desc.str(), inDatabank->GetVersion(), seqlength,
		chains.size(), kchain, usedChains, hits, res, outHSSP);
}

void CreateHSSP(
	CDatabankPtr		inDatabank,
	const MProtein&		inProtein,
	const fs::path&		inDataDir,
	const fs::path&		inFastaDir,
	const fs::path&		inJackHmmer,
	uint32				inIterations,
	vector<string>		inStockholmIds,
	ostream&			outHSSP)
{
	uint32 seqlength = 0;

	vector<hit_ptr> hits;
	vector<res_ptr> res;

	vector<mseq> alignments(inStockholmIds.size());
	vector<const MChain*> chains;

	uint32 kchain = 0;
	foreach (string ch, inStockholmIds)
	{
		if (ch.length() < 3 or ch[1] != '=')
			THROW(("Invalid chain/stockholm pair specified: '%s'", ch.c_str()));

		fs::path sfp = inDataDir / (ch.substr(2) + ".sto.bz2");
		if (not fs::exists(sfp))
		{
			// Stockholm file does not exist, create it
			const MChain& chain = inProtein.GetChain(ch[0]);
			string seq;
			chain.GetSequence(seq);

			RunJackHmmer(seq, inIterations, inFastaDir, inJackHmmer, inDatabank->GetID(), sfp);
		}

		fs::ifstream sf(sfp, ios::binary);
		if (not sf.is_open())
			THROW(("Could not open stockholm file '%s'", sfp.string().c_str()));

		io::filtering_stream<io::input> in;
		in.push(io::bzip2_decompressor());
		in.push(sf);

		ReadStockholm(in, alignments[kchain]);

		const MChain& chain = inProtein.GetChain(ch[0]);
		chains.push_back(&chain);
		
		// check to see if we need to 'cut' the alignment a bit
		// can happen if the stockholm file was created using a query
		// sequence that was a few residues longer than this chain.

		CheckAlignmentForChain(alignments[kchain], chains[kchain]);

		if (not res.empty())
			res.push_back(res_ptr(new ResidueHInfo(res.size() + 1)));

		ChainToHits(inDatabank, alignments[kchain], chain, hits, res);
		++kchain;
	}

	sort(hits.begin(), hits.end(), compare_hit());
	if (hits.size() > 9999)
		hits.erase(hits.begin() + 9999, hits.end());
	
	uint32 nr = 1;
	foreach (hit_ptr h, hits)
		h->nr = nr++;
	
	string usedChains;
	foreach (const MChain* chain, chains)
	{
		if (not usedChains.empty())
			usedChains += ',';
		usedChains += chain->GetChainID();
	}
	
	stringstream desc;
	desc	
	   << "HEADER     " + inProtein.GetHeader().substr(10, 40) << endl
	   << "COMPND     " + inProtein.GetCompound().substr(10) << endl
	   << "SOURCE     " + inProtein.GetSource().substr(10) << endl
	   << "AUTHOR     " + inProtein.GetAuthor().substr(10) << endl;

	CreateHSSPOutput(inProtein.GetID(), desc.str(), inDatabank->GetVersion(), seqlength,
		inProtein.GetChains().size(), kchain, usedChains, hits, res, outHSSP);
}

}
