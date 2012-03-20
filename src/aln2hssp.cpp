//  Copyright Maarten L. Hekkelman, Radboud University 2008.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include "mas.h"

#include <iostream>
#include <set>

#if P_UNIX
#include <wait.h>
#elif P_WIN
#include <conio.h>
#endif

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
#include <boost/program_options/config.hpp>
#include <boost/ptr_container/ptr_vector.hpp>
#include <boost/date_time/gregorian/gregorian.hpp>
#include <boost/date_time/date_clock_device.hpp>

#include "matrix.h"
#include "hssp.h"
#include "structure.h"
#include "utils.h"

using namespace std;
namespace ba = boost::algorithm;
namespace io = boost::iostreams;
namespace po = boost::program_options;
namespace fs = boost::filesystem;

// Globals section

fs::path gTempDir	= "/tmp/hssp-2/";
uint32 gMaxRunTime	= 3600;
uint32 gNrOfThreads;
int VERBOSE = 0;

#ifdef NDEBUG
int VERBOSE;
#endif

// main

int main(int argc, char* argv[])
{
#if P_UNIX
	// enable the dumping of cores to enable postmortem debugging
	rlimit l;
	if (getrlimit(RLIMIT_CORE, &l) == 0)
	{
		l.rlim_cur = l.rlim_max;
		if (l.rlim_cur == 0 or setrlimit(RLIMIT_CORE, &l) < 0)
			cerr << "Failed to set rlimit" << endl;
	}
#endif

	try
	{
		po::options_description desc("aln2hssp options");
		desc.add_options()
			("help,h",							 "Display help message")
			("input,i",		po::value<string>(), "Input PDB file (or PDB ID)")
			("chain",		po::value<vector<string>>(),
												 "Chain mapping in the form A=file-A.aln")
			("output,o",	po::value<string>(), "Output file, use 'stdout' to output to screen")
			("databank,b",	po::value<string>(), "Databank to use (default is uniprot)")
			("max-hits,m",	po::value<uint32>(), "Maximum number of hits to include (default = 1500)")
			("threshold",	po::value<float>(),  "Homology threshold adjustment (default = 0.05)")

			("verbose,v",						 "Verbose output")
			("debug,d",		po::value<int>(),	 "Debug level (for even more verbose output)")
			;
	
		po::positional_options_description p;
		p.add("input", 1);
		p.add("output", 2);
	
		po::variables_map vm;
		po::store(po::command_line_parser(argc, argv).options(desc).positional(p).run(), vm);

		//fs::path home = get_home();
		//if (fs::exists(home / ".mkhssprc"))
		//{
		//	fs::ifstream rc(home / ".mkhssprc");
		//	po::store(po::parse_config_file(rc, desc), vm);
		//}

		po::notify(vm);

		vector<string> mapped = vm["chain"].as<vector<string>>();
		if (vm.count("help") or not vm.count("input") or mapped.empty())
		{
			cerr << desc << endl;
			exit(1);
		}

		VERBOSE = vm.count("verbose") ? 1 : 0;
		if (vm.count("debug"))
			VERBOSE = vm["debug"].as<int>();
		
		string databank = "uniprot";
		if (vm.count("databank"))
			databank = vm["databank"].as<string>();
			
		uint32 maxhits = 1500;
		if (vm.count("max-hits"))
			maxhits= vm["max-hits"].as<uint32>();

		float threshold = 0.05f;
		if (vm.count("threshold"))
			threshold = vm["threshold"].as<float>();

		gNrOfThreads = boost::thread::hardware_concurrency();
//		if (vm.count("threads"))
//			gNrOfThreads = vm["threads"].as<uint32>();
//		if (gNrOfThreads < 1)
//			gNrOfThreads = 1;

		// what input to use
		string input = vm["input"].as<string>();
		io::filtering_stream<io::input> in;

		ifstream infile(input.c_str(), ios_base::in | ios_base::binary);
		istringstream indata;

		if (ba::ends_with(input, ".bz2"))
		{
			in.push(io::bzip2_decompressor());
			input.erase(input.length() - 4, string::npos);
		}
		else if (ba::ends_with(input, ".gz"))
		{
			in.push(io::gzip_decompressor());
			input.erase(input.length() - 3, string::npos);
		}
		in.push(infile);
		
		MProtein a(in);
		
		// then calculate the secondary structure
		a.CalculateSecondaryStructure();

		// Where to write our HSSP file to:
		// either to cout or an (optionally compressed) file.
		ofstream outfile;
		io::filtering_stream<io::output> out;

		if (vm.count("output") and vm["output"].as<string>() != "stdout")
		{
			string output = vm["output"].as<string>();
			outfile.open(output.c_str(), ios_base::out|ios_base::trunc|ios_base::binary);
			
			if (not outfile.is_open())
				throw runtime_error("could not create output file");
			
			if (ba::ends_with(output, ".bz2"))
				out.push(io::bzip2_compressor());
			else if (ba::ends_with(output, ".gz"))
				out.push(io::gzip_compressor());
			out.push(outfile);
		}
		else
			out.push(cout);

		CreateHSSPForAlignments(a, maxhits, threshold, mapped, out);
	}
	catch (exception& e)
	{
		cerr << e.what() << endl;
		exit(1);
	}

#if P_WIN && P_DEBUG
	cerr << "Press any key to quit application ";
	char ch = _getch();
	cerr << endl;
#endif
	
	return 0;
}
