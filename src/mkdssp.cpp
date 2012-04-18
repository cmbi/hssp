// Copyright Maarten L. Hekkelman, Radboud University 2008-2011.
//   Distributed under the Boost Software License, Version 1.0.
//       (See accompanying file LICENSE_1_0.txt or copy at    
//             http://www.boost.org/LICENSE_1_0.txt)      
//
// A DSSP reimplementation

#include "mas.h"

#if defined(_MSC_VER)
#include <conio.h>
#include <ctype.h>
#endif

#include <fstream>

#include <boost/program_options.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/copy.hpp>
#if defined USE_COMPRESSION
#include <boost/iostreams/filter/bzip2.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#endif
#include <boost/algorithm/string.hpp>

#include "dssp.h"
#include "structure.h"

using namespace std;
namespace po = boost::program_options;
namespace io = boost::iostreams;
namespace ba = boost::algorithm;

int MULTI_THREADED = 1, VERBOSE;

int main(int argc, char* argv[])
{
	try
	{
		po::options_description desc("DSSP " VERSION " options");
		desc.add_options()
			("help,h",							 "Display help message")
			("input,i",		po::value<string>(), "Input file")
			("output,o",	po::value<string>(), "Output file, use 'stdout' to output to screen")
			("verbose,v",						 "Verbose output")
			("version",							 "Print version")
			("debug,d",		po::value<int>(),	 "Debug level (for even more verbose output)")
			;
	
		po::positional_options_description p;
		p.add("input", 1);
		p.add("output", 2);
	
		po::variables_map vm;
		po::store(po::command_line_parser(argc, argv).options(desc).positional(p).run(), vm);
		po::notify(vm);

		if (vm.count("version"))
		{
			cout << "mkdssp version " VERSION << endl;
			exit(0);
		}

		if (vm.count("help") or not vm.count("input"))
		{
			cerr << desc << endl
				 << endl
				 << "Examples: " << endl
				 << endl
				 << "To calculate the secondary structure for the file 1crn.pdb and" << endl
				 << "write the result to a file called 1crn.dssp, you type:" << endl
				 << endl
				 << "  dssp.exe -i 1crn.pdb -o 1crn.dssp" << endl
				 << endl;
#if defined(_MSC_VER)
			cerr << endl
				 << "DSSP is a command line application, use the 'Command prompt' application" << endl
				 << "to start dssp.exe. You can find the 'Command prompt' in the Start menu:" << endl
				 << endl
				 << "Start => Accessories => Command prompt" << endl
				 << endl
				 << endl
				 << "Press any key to continue..." << endl;
			char ch = _getch();
#endif
			exit(1);
		}

		VERBOSE = vm.count("verbose") != 0;
		if (vm.count("debug"))
			VERBOSE = vm["debug"].as<int>();
		
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
	
		// and finally report the secondary structure in the DSSP format
		// either to cout or an (optionally compressed) file.
		if (vm.count("output"))
		{
			string output = vm["output"].as<string>();
			
			ofstream outfile(output.c_str(), ios_base::out|ios_base::trunc|ios_base::binary);
			if (not outfile.is_open())
				throw runtime_error("could not create output file");
			
			io::filtering_stream<io::output> out;
#if defined USE_COMPRESSION
			if (ba::ends_with(output, ".bz2"))
				out.push(io::bzip2_compressor());
			else if (ba::ends_with(output, ".gz"))
				out.push(io::gzip_compressor());
#endif
			out.push(outfile);
			
			WriteDSSP(a, out);
		}
		else
			WriteDSSP(a, cout);
	}
	catch (const exception& e)
	{
		cerr << "DSSP could not be created due to an error:" << endl
			 << e.what() << endl;
		exit(1);
	}
	
	return 0;
}

