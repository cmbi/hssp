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
#include <boost/ptr_container/ptr_vector.hpp>

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

extern "C" {
#include <clustal-omega.h>
}

using namespace std;
namespace ba = boost::algorithm;
namespace io = boost::iostreams;
namespace po = boost::program_options;

int nrOfThreads = boost::thread::hardware_concurrency();

struct insertion
{
	uint32		ipos, jpos;
	string		seq;
};

struct hit
{
				hit(const string& id, const string& seq) : nr(0), id(id), seq(seq) {}

	uint32		nr;
	string		id, acc, desc, pdb;
	string		seq;
	uint32		ifir, ilas, jfir, jlas, lali, ngap, lgap, lseq2;
	float		ide, wsim;
	uint32		identical, similar;
	vector<insertion>
				insertions;

	bool		operator<(const hit& rhs) const 	{ return ide > rhs.ide; }
};

ostream& operator<<(ostream& os, const hit& h)
{
	static boost::format fmt("%5.5d : %12.12s%4.4s    %4.2f  %4.2f %4.4d %4.4d %4.4d %4.4d %4.4d %4.4d %4.4d %4.4d ");
	
	os << fmt % h.nr % h.id % h.pdb
			  % h.ide % h.wsim % h.ifir % h.ilas % h.jfir % h.jlas % h.lali
			  % h.ngap % h.lgap % h.lseq2;
	
	return os;
}

typedef unique_ptr<hit> hit_ptr;

hit_ptr CreateHit(CDatabankPtr db, const string& id, const string& q, const string& s)
{
	hit_ptr result;
	
	assert(q.length() == s.length());

	sequence sq = encode(q);
	sequence ss = encode(s);

	const substitution_matrix m("BLOSUM62");
	
	uint32 b = 0, e = q.length();
	while (b < e)
	{
		if (q[b] == '-' or s[b] == '-' or m(sq[b], ss[b]) < 0)
		{
			++b;
			continue;
		}
		break;
	}

	while (b < e)
	{
		if (q[e - 1] == '-' or s[e - 1] == '-' or m(sq[e - 1], ss[e - 1]) <= 0)
		{
			--e;
			continue;
		}
		break;
	}
	
	result.reset(new hit(id, s));
	hit& h = *result;

	h.ifir = b + 1;
	for (uint32 i = 0; i < b and q[i] == '-'; ++i)
		--h.ifir;
	
	h.ilas = h.ifir + (e - b) - 1;
	for (uint32 i = e - 1; i > b and q[i] == '-'; --i)
		--h.ilas;

	h.jfir = b + 1;
	for (uint32 i = 0; i < b and s[i] == '-'; ++i)
		--h.jfir;

	h.jlas = h.jfir + (e - b) - 1;
	for (uint32 i = e - 1; i > b and s[i] == '-'; --i)
		--h.jlas;

	h.lseq2 = 0;
	h.lgap = 0;
	h.ngap = 0;
	h.identical = 0;
	h.similar = 0;
	
	bool gap = true;

	uint32 sb = 0;
	while (sb < s.length() and s[sb] == '-')
		++sb;
	
	uint32 se = s.length();
	while (se > sb and s[se - 1] == '-')
		--se;
	
	for (uint32 i = sb; i < se; ++i)
	{
		if (s[i] == '-' and q[i] == '-')
			continue;
		
		if (s[i] == '-')
		{
			if (not gap)
				++h.ngap;
			gap = true;
			++h.lgap;
		}
		else
		{
			++h.lseq2;
			gap = false;
		}
	}
	
	h.lali = 0;
	
	for (uint32 i = b; i < e; ++i)
	{
		if (q[i] == '-' and s[i] == '-')
		{
			--h.ilas;
			--h.jlas;
			continue;
		}

		++h.lali;

		if (q[i] == s[i])
		{
			++h.identical;
			++h.similar;
			continue;
		}
		
		if (q[i] == '-' or s[i] == '-')
		{
			if (s[i] == '-')
				--h.jlas;
			else
				--h.ilas;
			continue;
		}
		
		if (m(sq[i], ss[i]) > 0)
			++h.similar;
	}
	
	h.ide = float(h.identical) / float(h.lali);
	h.wsim = float(h.similar) / float(h.lali);
	
	return result;
}

void CreateHSSP(CDatabankPtr inDatabank, MProtein& inProtein,
	float identity_cutoff, opts_t& coo, ostream& os)
{
	stringstream dssp;
	WriteDSSP(inProtein, dssp);

	vector<string> seqs;
	inProtein.GetSequences(back_inserter(seqs));
	seqs.erase(unique(seqs.begin(), seqs.end()), seqs.end());

	// blast parameters
	float expect = 1.0;
	bool filter = true, gapped = true;
	int wordsize = 3, gapOpen = 11, gapExtend = 1, maxhits = 1500;
	string matrix = "BLOSUM62";
	
	foreach (const string& seq, seqs)
	{
		vector<uint32> hits;

		CDbAllDocIterator data(inDatabank.get());
		CBlast blast(seq, matrix, wordsize, expect, filter, gapped, gapOpen, gapExtend, maxhits);
		
		if (blast.Find(*inDatabank, data, nrOfThreads))
		{
			CBlastHitList hits(blast.Hits());
			
			// now create the alignment using clustalo
			
			mseq_t* msa;
			NewMSeq(&msa);
			vector<string> ids;
			
			AddSeq(&msa, const_cast<char*>(inProtein.GetID().c_str()),
				const_cast<char*>(seq.c_str()));
			
			foreach (const CBlastHit& hit, hits)
			{
				string seq, id;

				inDatabank->GetSequence(hit.DocumentNr(),
					inDatabank->GetSequenceNr(hit.DocumentNr(), hit.SequenceID()),
					seq);
				id = hit.DocumentID();
				
				ids.push_back(id);
				AddSeq(&msa, const_cast<char*>(id.c_str()), const_cast<char*>(seq.c_str()));
			}
			
			msa->seqtype = SEQTYPE_PROTEIN;
			msa->aligned = false;
			
			if (Align(msa, nil, &coo))
			{
				FreeMSeq(&msa);
				throw mas_exception("Fatal error creating alignment");
			}
			
			boost::ptr_vector<hit> hssp;
			for (int i = 1; i < msa->nseqs; ++i)
			{
				hit_ptr hit = CreateHit(inDatabank, msa->sqinfo[i].name,
					msa->seq[0], msa->seq[i]);
				if (hit->ide >= identity_cutoff)
					hssp.push_back(hit.release());
			}
			
			if (hssp.size() < msa->nseqs)	// repeat alignment with the new, smaller set of remaining hits
			{
				mseq_t* rs;
				NewMSeq(&rs);
				
				string s = seq;
				ba::erase_all(s, "-");
				
				AddSeq(&rs, const_cast<char*>(inProtein.GetID().c_str()),
					const_cast<char*>(s.c_str()));

				foreach (hit& h, hssp)
				{
					s = h.seq;
					ba::erase_all(s, "-");
					
					AddSeq(&rs, const_cast<char*>(h.id.c_str()), const_cast<char*>(s.c_str()));
				}

				rs->seqtype = SEQTYPE_PROTEIN;
				rs->aligned = false;
				
				if (Align(rs, nil, &coo))
				{
					FreeMSeq(&msa);
					FreeMSeq(&rs);
					throw mas_exception("Fatal error creating alignment");
				}
				
				FreeMSeq(&msa);
				msa = rs;
				
				hssp.clear();
				for (int i = 1; i < msa->nseqs; ++i)
				{
					hit_ptr hit = CreateHit(inDatabank, msa->sqinfo[i].name,
						msa->seq[0], msa->seq[i]);
					if (hit->ide >= identity_cutoff)
						hssp.push_back(hit.release());
				}
			}

			sort(hssp.begin(), hssp.end());
			uint32 nr = 1;
			foreach (hit& h, hssp)
			{
				h.nr = nr++;
				cout << h << endl;
			}
			
			cout << endl;
			foreach (hit& h, hssp)
				cout << h.seq << endl;

//			cout << endl;
//			for (int i = 0; i < msa->nseqs; ++i)
//			{
//				string id = msa->sqinfo[i].name;
//				id.insert(id.end(), 12 - id.length(), ' ');
//				
//				cout << id << ' ' << msa->seq[i] << endl;
//			}

			FreeMSeq(&msa);
		}
	}
	
	// finally create a HSSP file
	
	string hssp;
	
	
	
}

int main(int argc, char* argv[])
{
	try
	{
	    LogDefaultSetup(&rLog);

		opts_t coo;
		SetDefaultAlnOpts(&coo);
		
		po::options_description desc("DSSP options");
		desc.add_options()
			("help,h",							 "Display help message")
			("input,i",		po::value<string>(), "Input PDB file")
			("output,o",	po::value<string>(), "Output file, use 'stdout' to output to screen")
			("blastdb,b",	po::value<string>(), "Blast databank to use (default is uniprot)")
			("cutoff,c",	po::value<float>(),  "Identity cut off")
//			("maxhom",		po::value<string>(), "Path to the maxhom application")
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
		{
		    rLog.iLogLevelEnabled = LOG_DEBUG;
			VERBOSE = vm["debug"].as<int>();
		}
		
		string databank = "uniprot";
		if (vm.count("blastdb"))
			databank = vm["blastdb"].as<string>();
			
//		string maxhom = "maxhom";
//		if (vm.count("maxhom"))
//			maxhom = vm["maxhom"].as<string>();
		
		if (vm.count("threads"))
			nrOfThreads = vm["threads"].as<int>();
			
		float identity_cutoff = 0.4;
		if (vm.count("cutoff"))
			identity_cutoff = vm["cutoff"].as<float>();

		// init clustalo
		
		InitClustalOmega(nrOfThreads);

		CDatabankTable sDBTable;
		CDatabankPtr db = sDBTable.Load(databank);

//hit_ptr hit = CreateHit(db, "THNA_PHOLI",
//	"----------------------------TTCCPSIVARSNFNVCRLPGT-PEAICATYTGCIIIPGATCPGDYAN-------------------------",
//	"----------------------------KSCCPSTTARNIYNTCRLTGT-SRPTCASLSGCKIISGSTCBSGWBH-------------------------");
//		cout << *hit << endl;
//		return 0;
//
//
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
		CreateHSSP(db, a, identity_cutoff, coo, os);
		
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
