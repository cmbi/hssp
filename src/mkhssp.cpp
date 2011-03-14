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