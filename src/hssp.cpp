// A DSSP reimplementation
//
//	Copyright, M.L. Hekkelman, UMC St. Radboud, Nijmegen
//

#include "MRS.h"
#include "CSequence.h"

#include "mas.h"

#include <iostream>

#include <boost/program_options.hpp>
#include <boost/filesystem/fstream.hpp>
#include <boost/foreach.hpp>
#define foreach BOOST_FOREACH
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/copy.hpp>
#include <boost/iostreams/filter/bzip2.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/date_time/gregorian/gregorian.hpp>
#include <boost/date_time/date_clock_device.hpp>

#include "structure.h"
#include "ioseq.h"
#include "matrix.h"
#include "utils.h"

#include "CDatabank.h"
#include "CDatabankTable.h"
#include "CBlast.h"

using namespace std;
namespace fs = boost::filesystem;
namespace po = boost::program_options;
namespace io = boost::iostreams;
namespace ba = boost::algorithm;

int MULTI_THREADED = 1;

void alignForHSSP(vector<entry>& data, vector<entry*>& alignment)
{
	string matrix;
	float gop, gep, magic;
	
	matrix = "BLOSUM";
	gop = 10;
	gep = 0.2;
	magic = 0.1;

	substitution_matrix_family mat(matrix);
	joined_node node;
	
	progress pr("creating alignment", data.size() - 1);

	for (uint32 i = 1; i < data.size(); ++i)
	{
		entry ea(data.front());
		entry eb(data[i]);
		
		vector<entry*> a, b, c;
		a.push_back(&ea);
		b.push_back(&eb);
		
		align(&node, a, b, c, mat, gop, gep, magic, true);
		
		data[i].m_seq.clear();
		
		for (uint32 j = 0; j < ea.m_seq.length(); ++j)
		{
			if (ea.m_seq[j] != kSignalGapCode)
				data[i].m_seq += eb.m_seq[j];
		}
		
//		alignment.push_back(c.back());
		
		pr.step(1);
	}

	foreach (entry& e, data)
		alignment.push_back(&e);

//	alignment.push_back(&data.front());
//	for (uint32 i = 1; i < data.size(); ++i)
//	{
//		vector<entry*> b;
//		b.push_back(&data[i]);
//		
//		vector<entry*> c;
//		align(&node, alignment, b, c, mat, gop, gep, magic, true);
//		
//		alignment = c;
//		
//		pr.step(1);
//	}
}

void CreateHSSP(MProtein& protein, ostream& os)
{
	CDatabankTable dbTable;
	
	CDatabankPtr db = dbTable.Load("sprot");
	
	foreach (const MChain* chain, protein.GetChains())
	{
		char chainID = chain->GetChainID();

		vector<entry> data;
		entry e(data.size(), protein.GetID() + chainID);
		
		protein.GetSequence(chainID, e);
		
		data.push_back(e);
		
		CBlast blast(decode(e.m_seq), "BLOSUM62", 3, 10, true, true, 11, 1, 1500);
		
		CDbAllDocIterator all(db.get());
		if (blast.Find(*db, all, boost::thread::hardware_concurrency()))
		{
			CBlastHitIterator hitIter = blast.Hits();
			
			vector<uint32> hits;
			for (uint32 n = 1; hitIter.Next(); ++n)
				hits.push_back(hitIter.DocumentNr());
			
			if (VERBOSE)
				cerr << "blast resulted in " << hits.size() << " hits" << endl;
			
			foreach (uint32 hit, hits)
			{
				string id = db->GetDocumentID(hit);
//				cerr << id << endl;

				string seq;
				db->GetSequence(hit, 0, seq);
				entry e(data.size(), id, encode(seq));
				
				data.push_back(e);
			}
			
			vector<entry*> alignment;
			alignForHSSP(data, alignment);
			report(alignment, cout, "clustalw");
		}
	}
}

int main(int argc, char* argv[])
{
	try
	{
		po::options_description desc("DSSP options");
		desc.add_options()
			("help,h",							 "Display help message")
			("input,i",		po::value<string>(), "Input file")
			("output,o",	po::value<string>(), "Output file, use 'stdout' to output to screen")
			("verbose,v",						 "Verbose output")
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
	
		// and finally report the secondary structure in the DSSP format
		// either to cout or an (optionally compressed) file.
		if (vm.count("output"))
		{
			string output = vm["output"].as<string>();
			
			ofstream outfile(output.c_str(), ios_base::out|ios_base::trunc|ios_base::binary);
			if (not outfile.is_open())
				throw runtime_error("could not create output file");
			
			io::filtering_stream<io::output> out;
			if (ba::ends_with(output, ".bz2"))
				out.push(io::bzip2_compressor());
			else if (ba::ends_with(output, ".gz"))
				out.push(io::gzip_compressor());
			
			out.push(outfile);
			
			CreateHSSP(a, out);
		}
		else
			CreateHSSP(a, cout);
	}
	catch (const exception& e)
	{
		cerr << e.what() << endl;
		exit(1);
	}
	
	return 0;
}

