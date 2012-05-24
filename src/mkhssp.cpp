//  Copyright Maarten L. Hekkelman, Radboud University 2008.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include "mas.h"

#if defined(_MSC_VER)
#include <conio.h>
#endif

#include <boost/algorithm/string.hpp>
#include <boost/filesystem/fstream.hpp>
#include <boost/program_options.hpp>
#include <boost/program_options/config.hpp>
#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/filter/zlib.hpp>
#include <boost/iostreams/device/back_inserter.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/copy.hpp>
#include <boost/iostreams/filter/bzip2.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/foreach.hpp>
#define foreach BOOST_FOREACH

#include "utils.h"
#include "structure.h"
#include "hssp-nt.h"

using namespace std;
namespace ba = boost::algorithm;
namespace fs = boost::filesystem;
namespace io = boost::iostreams;
namespace po = boost::program_options;

// --------------------------------------------------------------------

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
		po::options_description desc("MKHSSP options");
		desc.add_options()
			("help,h",							 "Display help message")
			("input,i",		po::value<string>(), "Input PDB file (or PDB ID)")
			("output,o",	po::value<string>(), "Output file, use 'stdout' to output to screen")
			("databank,d",	po::value<vector<string>>(),
												 "Databank(s) to use")
			("threads,a",	po::value<uint32>(), "Number of threads (default is maximum)")
			("use-seqres",	po::value<bool>(),	 "Use SEQRES chain instead of chain based on ATOM records (values are true of false, default is true)")
			("min-length",	po::value<uint32>(), "Minimal chain length")
			("fragment-cutoff",
							po::value<float>(),  "Minimal alignment length as fraction of chain length (default = 0.75)")
			("gap-open,O",	po::value<float>(),  "Gap opening penalty (default is 30.0)")
			("gap-extend,E",po::value<float>(),  "Gap extension penalty (default is 2.0)")
			("threshold",	po::value<float>(),  "Homology threshold adjustment (default = 0.05)")
			("max-hits,m",	po::value<uint32>(), "Maximum number of hits to include (default = 1500)")
			("fetch-dbrefs",					 "Fetch DBREF records for each UniProt ID")
			("verbose,v",						 "Verbose output")
			;
	
		po::positional_options_description p;
		p.add("input", 1);
		p.add("output", 2);
	
		po::variables_map vm;
		po::store(po::command_line_parser(argc, argv).options(desc).positional(p).run(), vm);

		fs::path home = get_home();
		if (fs::exists(home / ".mkhssprc"))
		{
			fs::ifstream rc(home / ".mkhssprc");
			po::store(po::parse_config_file(rc, desc), vm);
		}

		po::notify(vm);

		if (vm.count("help") or not vm.count("input") or vm.count("databank") == 0)
		{
			cerr << desc << endl;
			exit(1);
		}

		VERBOSE = vm.count("verbose") > 0;
		
		vector<fs::path> databanks;
		vector<string> dbs = vm["databank"].as<vector<string>>(); 
		foreach (string db, dbs)
		{
			databanks.push_back(db);
			if (not fs::exists(databanks.back()))
				throw mas_exception(boost::format("Databank %s does not exist") % db);
		}
		
		bool useSeqRes = true;
		if (vm.count("use-seqres"))
			useSeqRes = vm["use-seqres"].as<bool>();
		
		uint32 minlength = 25;
		if (vm.count("min-length"))
			minlength= vm["min-length"].as<uint32>();

		uint32 maxhits = 5000;
		if (vm.count("max-hits"))
			maxhits= vm["max-hits"].as<uint32>();

		float gapOpen = 30;
		if (vm.count("gap-open"))
			gapOpen = vm["gap-open"].as<float>();

		float gapExtend = 2;
		if (vm.count("gap-extend"))
			gapExtend = vm["gap-extend"].as<float>();

		float threshold = HSSP::kThreshold;
		if (vm.count("threshold"))
			threshold = vm["threshold"].as<float>();

		float fragmentCutOff = HSSP::kFragmentCutOff;
		if (vm.count("fragment-cutoff"))
			fragmentCutOff = vm["fragment-cutoff"].as<float>();

		bool fetchDbRefs = vm.count("fetch-dbrefs") > 0;

		uint32 threads = boost::thread::hardware_concurrency();
		if (vm.count("threads"))
			threads = vm["threads"].as<uint32>();
		if (threads < 1)
			threads = 1;
			
		// what input to use
		string input = vm["input"].as<string>();
		io::filtering_stream<io::input> in;
		ifstream infile(input.c_str(), ios_base::in | ios_base::binary);
		if (not infile.is_open())
			throw runtime_error("Error opening input file");

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

		// read protein and calculate the secondary structure
		MProtein a(in);
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

		// create the HSSP file
		HSSP::CreateHSSP(a, databanks, maxhits, minlength,
			gapOpen, gapExtend, threshold, fragmentCutOff, threads, fetchDbRefs, out);
	}
	catch (exception& e)
	{
		cerr << e.what() << endl;
		exit(1);
	}

//#if defined(_MSC_VER) && ! NDEBUG
//	cerr << "Press any key to quit application ";
//	char ch = _getch();
//	cerr << endl;
//#endif
	
	return 0;
}

