//  Copyright Maarten L. Hekkelman, Radboud University 2008.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include "MRS.h"

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
#include <boost/ptr_container/ptr_vector.hpp>
#include <boost/date_time/gregorian/gregorian.hpp>
#include <boost/date_time/date_clock_device.hpp>

#include "CDatabank.h"
#include "CDatabankTable.h"
#include "CBlast.h"
#include "CQuery.h"

#include "mas.h"
#include "matrix.h"
#include "dssp.h"
#include "blast.h"
#include "structure.h"
#include "utils.h"
#include "hmmer-hssp.h"

using namespace std;
namespace ba = boost::algorithm;
namespace io = boost::iostreams;
namespace po = boost::program_options;

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
			("input,i",		po::value<string>(), "Input PDB file")
			("output,o",	po::value<string>(), "Output file, use 'stdout' to output to screen")
			("databank,b",	po::value<string>(), "Databank to use (default is uniprot)")
			("fastadir,f",	po::value<string>(), "Directory containing fasta databank files)")
			("jackhmmer",	po::value<string>(), "Jackhmmer executable path (default=/usr/local/bin/jackhmmer)")
			("iterations",	po::value<uint32>(), "Number of jackhmmer iterations (default = 5)")
			("dssp-id",		po::value<string>(), "Calculate HSSP for a specific DSSP ID")

			("datadir",		po::value<string>(), "Data directory containing stockholm files")
			("chain",		po::value<vector<string>>(),
												 "Mappings for chain => stockholm file")
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
		if (vm.count("databank"))
			databank = vm["databank"].as<string>();
			
		vector<string> chains;
		if (vm.count("chain"))
			chains = vm["chain"].as<vector<string>>();

		fs::path jackhmmer("/usr/local/bin/jackhmmer");
		if (vm.count("jackhmmer"))
			jackhmmer = fs::path(vm["jackhmmer"].as<string>());
		if (chains.empty() and not fs::exists(jackhmmer))
			throw mas_exception("Jackhmmer executable not found");
			
		fs::path fastadir("/data/fasta");
		if (vm.count("fastadir"))
			fastadir = fs::path(vm["fastadir"].as<string>());
		if (chains.empty() and not fs::exists(fastadir))
			throw mas_exception("Fasta databank directory not found");
			
		uint32 iterations = 5;
		if (vm.count("iterations"))
			iterations = vm["iterations"].as<uint32>();

		fs::path datadir(".");
		if (vm.count("datadir"))
			datadir = fs::path(vm["datadir"].as<string>());
		if (not fs::exists(datadir))
			throw mas_exception("Data directory not found");

		// got parameters

		CDatabankTable sDBTable;
		CDatabankPtr db = sDBTable.Load(databank);

		// what input to use
		string input = vm["input"].as<string>();

		ifstream infile(input.c_str(), ios_base::in | ios_base::binary);
		if (not infile.is_open())
			throw runtime_error("No such file");
		
		io::filtering_stream<io::input> in;
		if (ba::ends_with(input, ".bz2"))
			in.push(io::bzip2_decompressor());
		else if (ba::ends_with(input, ".gz"))
			in.push(io::gzip_decompressor());
		in.push(infile);

		// OK, we've got the file, now create a protein
		MProtein a(in);
		
		// then calculate the secondary structure
		a.CalculateSecondaryStructure();

		// Where to write our HSSP file to:
		// either to cout or an (optionally compressed) file.
		if (vm.count("output") and vm["output"].as<string>() != "stdout")
		{
			string output = vm["output"].as<string>();
			
			ofstream outfile(output.c_str(), ios_base::out|ios_base::trunc|ios_base::binary);
			if (not outfile.is_open())
				throw runtime_error("could not create output file");
			
			try
			{
				io::filtering_stream<io::output> out;
				if (ba::ends_with(output, ".bz2"))
					out.push(io::bzip2_compressor());
				else if (ba::ends_with(output, ".gz"))
					out.push(io::gzip_compressor());
				out.push(outfile);
	
				// and the final HSSP file
				if (chains.empty())
					hmmer::CreateHSSP(db, a, fastadir, jackhmmer, iterations, 25, out);
				else
					hmmer::CreateHSSP(db, a, datadir, chains, out);
			}
			catch (...)
			{
				outfile.close();
				fs::remove(output);
				throw;
			}
		}
		else
		{
			if (chains.empty())
				hmmer::CreateHSSP(db, a, fastadir, jackhmmer, iterations, 25, cout);
			else
				hmmer::CreateHSSP(db, a, datadir, chains, cout);
		}
	}
	catch (exception& e)
	{
		cerr << e.what() << endl;
		exit(1);
	}

#if P_WIN && P_DEBUG
	cerr << "Press any key to quit application ";
	char ch = _getch();
#endif
	
	return 0;
}
