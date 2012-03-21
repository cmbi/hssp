//  Copyright Maarten L. Hekkelman, Radboud University 2008.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include "mas.h"
#include "MRS.h"

#if P_WIN
#include <conio.h>
#endif

#include <iostream>
#include <fstream>

#include <boost/program_options.hpp>

#include "hmmer-hssp.h"
#include "CDatabankTable.h"

using namespace std;
namespace po = boost::program_options;
namespace fs = boost::filesystem;

// Globals section

CDatabankTable gDBTable;
fs::path gTempDir	= "/tmp/hssp-2/";
uint32 gMaxRunTime	= 3600;
uint32 gNrOfThreads;

// main

int main(int argc, char* argv[])
{
	try
	{
		po::options_description desc("sto2fa options");
		desc.add_options()
			("help,h",							 "Display help message")
			//("query,q",		po::value<string>(), "Input FastA formatted query file")
			("input,i",		po::value<string>(), "Input Stockholm file")
			("output,o",	po::value<string>(), "Output FastA file")
			("verbose,v",						 "Verbose mode")
			;
	
		po::positional_options_description p;
		//p.add("query", 1);
		p.add("input", 1);
		p.add("output", 2);
	
		po::variables_map vm;
		po::store(po::command_line_parser(argc, argv).options(desc).positional(p).run(), vm);
		po::notify(vm);

		if (vm.count("help") or not vm.count("input") or not vm.count("output"))
		{
			cerr << desc << endl;
			exit(1);
		}

		VERBOSE = vm.count("verbose") ? 1 : 0;
		if (vm.count("debug"))
			VERBOSE = vm["debug"].as<int>();
		
//		uint32 maxhits = 1500;
//		if (vm.count("max-hits"))
//			maxhits= vm["max-hits"].as<uint32>();
//

		// got parameters
		// what input to use
		
		//string queryFile = vm["query"].as<string>();
		//
		//ifstream f(queryFile);
		//if (not f.is_open())
		//	THROW(("Could not open query file"));
		//
		//string line;
		//getline(f, line);
		
		
		string input = vm["input"].as<string>();
		string output = vm["output"].as<string>();
		
		hmmer::ConvertHmmerAlignment("", input, output);
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
