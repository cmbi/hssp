//  Copyright Maarten L. Hekkelman, Radboud University 2008.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include "MRS.h"

#include <iostream>
#include <set>
#include <wait.h>

#include <boost/lexical_cast.hpp>
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
#include <boost/iostreams/filter/bzip2.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/foreach.hpp>
#define foreach BOOST_FOREACH
#include <boost/program_options.hpp>

#include "CDatabank.h"
#include "CDatabankTable.h"
#include "CBlast.h"
#include "CQuery.h"

#include "mas.h"
#include "dssp.h"
#include "structure.h"

using namespace std;
namespace ba = boost::algorithm;
namespace io = boost::iostreams;
namespace po = boost::program_options;

int BLAST_THREADS = 0;

void BlastSequence(
	CDatabankPtr		inDatabank,
	const string&		inSequence,
	vector<uint32>&		outHits)
{
	float expect = 1.0;
	bool filter = true, gapped = true;
	int wordsize = 3, gapOpen = 11, gapExtend = 1, threads, maxhits = 1500;
	string matrix = "BLOSUM62";
	
	threads = BLAST_THREADS;
	if (threads < 1)
		threads = boost::thread::hardware_concurrency();

	auto_ptr<CDocIterator> data(new CDbAllDocIterator(inDatabank.get()));
	CBlast blast(inSequence, matrix, wordsize, expect, filter, gapped, gapOpen, gapExtend, maxhits);
	
	if (blast.Find(*inDatabank, *data, threads))
	{
		CBlastHitIterator hits = blast.Hits();
		
		while (hits.Next())
			outHits.push_back(hits.DocumentNr());
	}
}

void WriteToFD(
	int						inFD,
	string					inText)
{
	inText += '\n';
	
	for (;;)
	{
		int r = write(inFD, inText.c_str(), inText.length());
		if (r == -1 and errno == EAGAIN)
			continue;
		if (r != static_cast<int>(inText.length()))
			THROW(("Failed to write maxhom command"));
		break;
	}		 
}

void GetHSSPForHitsAndDSSP(
	CDatabankPtr			inDatabank,
	const string&			inMaxHom,
	const string&			inPDBID,
	const vector<uint32>&	inHits,
	const string&			inDSSP,
	int						inMaxAlign,
	ostream&				outHSSP)
{
	int threshold = 5;
	
	HUuid uuid;
	
	fs::path rundir("/data/tmp/hssp-2/");
	rundir /= boost::lexical_cast<string>(uuid);
	fs::create_directories(rundir);
	
	// create 'hits' file containing swiss-prot entries for all hits
	fs::ofstream hits(rundir / "hits");
	for (vector<uint32>::const_iterator hit = inHits.begin(); hit != inHits.end(); ++hit)
		hits << inDatabank->GetDocument(*hit);
	hits.close();
	
	// write out the dssp file
	fs::ofstream dssp(rundir / (inPDBID + ".dssp"));
	dssp << inDSSP;
	dssp.close();
	
	// start a maxhom
	int ifd[2];
	
	if (pipe(ifd) != 0)
		THROW(("Failed to create pipe: %s", strerror(errno)));
	
	int pid = fork();
	
	if (pid == -1)
	{
		close(ifd[0]);
		close(ifd[1]);
		
		THROW(("fork failed: %s", strerror(errno)));
	}
	
	if (pid == 0)	// the child process (will be maxhom)
	{
		fs::current_path(rundir);
		
		setpgid(0, 0);
		
		dup2(ifd[0], STDIN_FILENO);
		close(ifd[0]);
		close(ifd[1]);
		
		int fd = open("maxhom.log", O_CREAT | O_RDWR, 0666);
		dup2(fd, STDOUT_FILENO);
		dup2(fd, STDERR_FILENO);
		
		char* argv[] = {
			const_cast<char*>(inMaxHom.c_str()),
			nil
		};
		
		(void)execve(inMaxHom.c_str(), argv, environ);
	}
	
	close(ifd[0]);
	int fd = ifd[1];
	
	// OK, so now we pass maxhom the commands and wait for it to exit
	WriteToFD(fd, "COMMAND NO");
	WriteToFD(fd, "BATCH");
	WriteToFD(fd, (boost::format("PID: %d") % getpid()).str());
	WriteToFD(fd, "SEQ_1 " + inPDBID + ".dssp");
	WriteToFD(fd, "SEQ_2 hits");
	WriteToFD(fd, "2_PROFILES NO");
	WriteToFD(fd, "METRIC LACHLAN");
	WriteToFD(fd, "NORM_PROFILE DISABLED");
	WriteToFD(fd, "MEAN_PROFILE ignored");
	WriteToFD(fd, "FACTOR_GAPS ignored");
	WriteToFD(fd, "SMIN -0.50");
	WriteToFD(fd, "SMAX 1.00");
	WriteToFD(fd, "GAP_OPEN 3.00");
	WriteToFD(fd, "GAP_ELONG 0.10");
	WriteToFD(fd, "WEIGHT1 YES");
	WriteToFD(fd, "WEIGHT2 NO");
	WriteToFD(fd, "WAY3-ALIGN NO");
	WriteToFD(fd, "INDEL_1 YES");
	WriteToFD(fd, "INDEL_2 YES");
	WriteToFD(fd, "RELIABILITY NO");
	WriteToFD(fd, "FILTER_RANGE 10.0");
	WriteToFD(fd, "NBEST 1");
	WriteToFD(fd, (boost::format("MAXALIGN %d") % inMaxAlign).str());
	WriteToFD(fd, (boost::format("THRESHOLD FORMULA%+d") % threshold).str());
	WriteToFD(fd, "SORT DISTANCE");
	WriteToFD(fd, "HSSP out.hssp");
	WriteToFD(fd, "SAME_SEQ_SHOW YES");
	WriteToFD(fd, "SUPERPOS NO");
	WriteToFD(fd, "PDB_PATH OPTION DISABLED");
	WriteToFD(fd, "PROFILE_OUT NO");
	WriteToFD(fd, "STRIP_OUT NO");
	WriteToFD(fd, "LONG_OUT NO");
	WriteToFD(fd, "DOT_PLOT NO");
	WriteToFD(fd, "RUN");
	
	int status;
	waitpid(pid, &status, 0);
	
	if (status != 0)
		THROW(("maxhom exited with status %d", status));
	
	if (not fs::exists(rundir / "out.hssp"))
		THROW(("Maxhom failed to create an alignment"));

	// OK, got it! Read in the result and exit
	fs::ifstream result(rundir / "out.hssp");
	
	// update the SEQBASE line while copying over the data
	string dbVersion = inDatabank->GetVersion();
	
	for (;;)
	{
		string line;
		getline(result, line);
		if (result.eof())
			break;
		
		if (ba::starts_with(line, "SEQBASE    "))
		{
			line.erase(11, string::npos);
			line += dbVersion;
		}
		
		outHSSP << line << endl;
	}
}

void CreateHSSP(
	CDatabankPtr		inDatabank,
	const string&		inMaxHom,
	MProtein&			inProtein,
	ostream&			outHSSP)
{
	stringstream dssp;
	WriteDSSP(inProtein, dssp);

	vector<string> seqs;
	inProtein.GetSequences(back_inserter(seqs));
	seqs.erase(unique(seqs.begin(), seqs.end()), seqs.end());
	
	vector<uint32> hits;
	foreach (const string& seq, seqs)
	{
		vector<uint32> h1;
		BlastSequence(inDatabank, seq, h1);
		sort(h1.begin(), h1.end());
		
		vector<uint32> h2;
		set_union(h1.begin(), h1.end(), hits.begin(), hits.end(), back_inserter(h2));
		swap(h2, hits);
	}
	
	if (hits.empty())
		THROW(("No blast hits found for %s", inProtein.GetID().c_str()));
	
	string hssp;
	GetHSSPForHitsAndDSSP(inDatabank, inMaxHom, inProtein.GetID(), hits, dssp.str(), 2000, outHSSP);
}

int main(int argc, char* argv[])
{
	try
	{
		po::options_description desc("DSSP options");
		desc.add_options()
			("help,h",							 "Display help message")
			("input,i",		po::value<string>(), "Input PDB file")
			("output,o",	po::value<string>(), "Output file, use 'stdout' to output to screen")
			("blastdb,b",	po::value<string>(), "Blast databank to use (default is uniprot)")
			("maxhom",		po::value<string>(), "Path to the maxhom application")
			("threads,a",	po::value<int>(),	 "Number of threads to use (default is nr of CPU's)")
			("verbose,v",						 "Verbose output")
			("debug,d",		po::value<int>(),	 "Debug level (for even more verbose output)")
			;
	
		po::positional_options_description p;
		p.add("input", 1);
		p.add("output", 2);
	
		po::variables_map vm;
		po::store(po::command_line_parser(argc, argv).options(desc).positional(p).run(), vm);
		po::notify(vm);

		if (vm.count("help") or not vm.count("input"))
		{
			cerr << desc << endl;
			exit(1);
		}

		VERBOSE = vm.count("verbose");
		if (vm.count("debug"))
			VERBOSE = vm["debug"].as<int>();
		
		string databank = "uniprot";
		if (vm.count("blastdb"))
			databank = vm["blastdb"].as<string>();
			
		string maxhom = "maxhom";
		if (vm.count("maxhom"))
			maxhom = vm["maxhom"].as<string>();
		
		if (vm.count("threads"))
			BLAST_THREADS = vm["threads"].as<int>();

		CDatabankTable sDBTable;
		CDatabankPtr db = sDBTable.Load(databank);

		// what input to use
		string input = vm["input"].as<string>();

		ifstream infile(input.c_str(), ios_base::in | ios_base::binary);
		if (not infile.is_open())
			throw runtime_error("No such file");
		
		io::filtering_stream<io::input> in;
		
#if defined USE_COMPRESSION
		if (ba::ends_with(input, ".bz2"))
			in.push(io::bzip2_decompressor());
		else if (ba::ends_with(input, ".gz"))
			in.push(io::gzip_decompressor());
#endif
		
		in.push(infile);

		// OK, we've got the file, now create a protein
		MProtein a(in);
		
		// then calculate the secondary structure
		a.CalculateSecondaryStructure();

		// and the final HSSP file
		// (we use a temporary stringstream, to avoid
		// creating empty files if something goes wrong.
		vector<char> hssp;
		
		io::filtering_ostream os(io::back_inserter(hssp));
		CreateHSSP(db, maxhom, a, os);
		
		io::filtering_istream is(boost::make_iterator_range(hssp));
		
		// Where to write our HSSP file to:
		// either to cout or an (optionally compressed) file.
		if (vm.count("output"))
		{
			string output = vm["output"].as<string>();
			
			ofstream outfile(output.c_str(), ios_base::out|ios_base::trunc|ios_base::binary);
			if (not outfile.is_open())
				throw runtime_error("could not create output file");
			
#if defined USE_COMPRESSION
			io::filtering_stream<io::output> out;
			if (ba::ends_with(output, ".bz2"))
				out.push(io::bzip2_compressor());
			else if (ba::ends_with(output, ".gz"))
				out.push(io::gzip_compressor());
#endif
			
			out.push(outfile);
			
			io::copy(is, out);
		}
		else
			io::copy(is, cout);
	}
	catch (exception& e)
	{
		cerr << e.what() << endl;
		exit(1);
	}
	
	return 0;
}
