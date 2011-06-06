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

extern "C" {
#include "clustal-omega.h"
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
	string		seq, saln;
	uint32		ifir, ilas, jfir, jlas, lali, ngap, lgap, lseq2;
	float		ide, wsim;
	uint32		identical, similar;
	vector<insertion>
				insertions;

	bool		operator<(const hit& rhs) const 	{ return ide > rhs.ide; }
	
	bool		IdentityAboveThreshold() const;
};

typedef shared_ptr<hit> hit_ptr;

struct ResidueHInfo
{
	const MResidue*	res;
	mseq_t*			msa;
	uint32			pos;
	uint32			seqNo, pdbNo;
	uint32			nocc, ndel, nins;
	float			entropy, weight;
	uint32			relent;
	uint32			var;
	uint32			dist[20];
};

typedef shared_ptr<ResidueHInfo>	res_ptr;

struct compare_hit
{
	bool operator()(hit_ptr a, hit_ptr b) const { return *a < *b; }
};

bool hit::IdentityAboveThreshold() const
{
	static vector<double> kHomologyThreshold;
	if (kHomologyThreshold.empty())
	{
		kHomologyThreshold.reserve(71);
		for (uint32 i = 10; i <= 80; ++i)
			kHomologyThreshold.push_back(2.9015 * pow(i, -0.562) + 0.05);
	}
	
	uint32 l = lali;
	if (l < 10)
		l = 0;
	else if (l > 80)
		l = 70;
	else
		l -= 10;
	
	assert(l < kHomologyThreshold.size());
	return kHomologyThreshold[l] < ide;
}

ostream& operator<<(ostream& os, const hit& h)
{
	static boost::format fmt("%5.5d : %12.12s%4.4s    %4.2f  %4.2f %4.4d %4.4d %4.4d %4.4d %4.4d %4.4d %4.4d %4.4d ");
	
	os << fmt % h.nr % h.id % h.pdb
			  % h.ide % h.wsim % h.ifir % h.ilas % h.jfir % h.jlas % h.lali
			  % h.ngap % h.lgap % h.lseq2;
	
	return os;
}

hit_ptr CreateHit(const string& id, const string& q, const string& s)
{
	hit_ptr result;
	
	assert(q.length() == s.length());

	sequence sq = encode(q);
	sequence ss = encode(s);

	// first remove common gaps
	sequence::iterator qi = sq.begin(), qj = qi, si = ss.begin(), sj = si;
	while (qi != sq.end())
	{
		if (*qi == kSignalGapCode and *si == kSignalGapCode)
		{
			++qi;
			++si;
			continue;
		}
		
		*qj++ = *qi++;
		*sj++ = *si++;
	}
	
	sq.erase(qj, sq.end());
	ss.erase(sj, ss.end());
	
//	if (VERBOSE)
//	{
//		cerr << "a: " << decode(sq) << endl
//			 << "b: " << decode(ss) << endl;
//	}
//
	const substitution_matrix m("BLOSUM62");
	
	sequence::iterator qb = sq.begin(), qe = sq.end(),
					   sb = ss.begin(), se = ss.end();

	result.reset(new hit(id, s));
	hit& h = *result;
	h.lseq2 = ss.length();
	h.lgap = h.ngap = h.identical = h.similar = 0;

	h.ifir = h.jfir = 1;
	h.ilas = h.jlas = sq.length();
	
	while (qb != qe)
	{
		if (*qb == kSignalGapCode)
		{
			++h.jfir;
			--h.ilas;
		}
		else if (*sb == kSignalGapCode)
		{
			++h.ifir;
			--h.jlas;
			--h.lseq2;
		}
		else if (m(*qb, *sb) <= 0)
		{
			++h.ifir;
			++h.jfir;
		}
		else
			break;
		
		++qb;
		++sb;
	}
	
	sq.erase(sq.begin(), qb); qb = sq.begin(); qe = sq.end();
	ss.erase(ss.begin(), sb); sb = ss.begin(); se = ss.end();

	while (qe != qb and (*(qe - 1) == kSignalGapCode or *(se - 1) == kSignalGapCode or m(*(qe - 1), *(se - 1)) <= 0))
	{
		if (*(se - 1) == kSignalGapCode)
			--h.lseq2;

		--h.ilas;
		--h.jlas;
		--qe;
		--se;
	}
	
	sq.erase(qe, sq.end());	qb = sq.begin(); qe = sq.end();
	ss.erase(se, ss.end());	sb = ss.begin(); se = ss.end();

//	if (VERBOSE)
//	{
//		cerr << "a: " << decode(sq) << endl
//			 << "b: " << decode(ss) << endl;
//	}
	
	h.lali = ss.length();
	h.saln = string(h.ifir - 1, ' ');
	
	bool gap = true;
	for (sequence::iterator si = sb, qi = qb; si != se; ++si, ++qi)
	{
		if (*si == kSignalGapCode)
		{
			if (not gap)
				++h.ngap;
			gap = true;
			++h.lgap;
			--h.lseq2;
			
			h.saln += ' ';
		}
		else if (*qi == kSignalGapCode)
		{
			if (not gap)
				++h.ngap;
			gap = true;
			++h.lgap;
		}
		else
		{
			gap = false;

			if (*qi == *si)
			{
				++h.identical;
				++h.similar;
			}
			else if (m(*qi, *si) > 0)
				++h.similar;

			h.saln += kAA[*si];
		}
	}

	h.ide = float(h.identical) / float(h.lali);
	h.wsim = float(h.similar) / float(h.lali);

//	if (VERBOSE)
//	{
//		cerr << h << endl;
//	}
	
	return result;
}

res_ptr CreateResidueHInfo(const MResidue* res, mseq_t* msa, uint32 pos)
{
	res_ptr r(new ResidueHInfo);
	
	r->res = res;
	r->msa = msa;
	r->pos = pos;
	
	r->nocc = r->ndel = r->nins = 0;
	fill(r->dist, r->dist + 20, 0);
	
	for (uint32 i = 0; i < msa->nseqs; ++i)
	{
		++r->nocc;
		// V   L   I   M   F   W   Y   G   A   P   S   T   C   H   R   K   Q   E   N   D
		switch (msa->seq[i][pos])
		{
			case 'V':	++r->dist[0]; break;
			case 'L':	++r->dist[1]; break;
			case 'I':	++r->dist[2]; break;
			case 'M':	++r->dist[3]; break;
			case 'F':	++r->dist[4]; break;
			case 'W':	++r->dist[5]; break;
			case 'Y':	++r->dist[6]; break;
			case 'G':	++r->dist[7]; break;
			case 'A':	++r->dist[8]; break;
			case 'P':	++r->dist[9]; break;
			case 'S':	++r->dist[10]; break;
			case 'T':	++r->dist[11]; break;
			case 'C':	++r->dist[12]; break;
			case 'H':	++r->dist[13]; break;
			case 'R':	++r->dist[14]; break;
			case 'K':	++r->dist[15]; break;
			case 'Q':	++r->dist[16]; break;
			case 'E':	++r->dist[17]; break;
			case 'N':	++r->dist[18]; break;
			case 'D':	++r->dist[19]; break;
			default:	--r->nocc; break;
		}
	}
	
	for (uint32 a = 0; a < 20; ++a)
		r->dist[a] = (100 * r->dist[a]) / r->nocc;
	
	
	//struct ResidueHInfo
	//{
	//	MResidue*		res;
	//	mseq_t*			msa;
	//	uint32			msa_pos;
	//	uint32			seqNo, pdbNo;
	//	uint32			nocc, ndel, nins;
	//	float			entropy, weight;
	//	uint32			relent;
	//	uint32			var;
	//	uint32			dist[20];
	//};

	return r;
}

void CreateHSSP(CDatabankPtr inDatabank, MProtein& inProtein,
	float identity_cutoff, opts_t& coo, ostream& os)
{
	stringstream dssp;
	WriteDSSP(inProtein, dssp);

#pragma warning("ketennamen opslaan")

//	seqs.erase(unique(seqs.begin(), seqs.end()), seqs.end());

	// blast parameters
	float expect = 1.0;
	bool filter = true, gapped = true;
	int wordsize = 3, gapOpen = 11, gapExtend = 1, maxhits = 1500;
	string matrix = "BLOSUM62";
	
	vector<hit_ptr> hssp;
	vector<res_ptr> result;
	
	foreach (const MChain* chain, inProtein.GetChains())
	{
		const vector<MResidue*>& residues(chain->GetResidues());
		
		string seq;
		chain->GetSequence(seq);
		
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
			
			vector<hit_ptr> c_hssp;
			for (int i = 1; i < msa->nseqs; ++i)
			{
				hit_ptr hit = CreateHit(msa->sqinfo[i].name, msa->seq[0], msa->seq[i]);

				if (hit->IdentityAboveThreshold())
					c_hssp.push_back(hit);
				else if (VERBOSE)
					cerr << "dropping " << hit->id << endl;
			}
			
			if (c_hssp.size() < msa->nseqs)	// repeat alignment with the new, smaller set of remaining hits
			{
				mseq_t* rs;
				NewMSeq(&rs);
				
				string s = seq;
				ba::erase_all(s, "-");
				
				AddSeq(&rs, const_cast<char*>(inProtein.GetID().c_str()),
					const_cast<char*>(s.c_str()));

				foreach (hit_ptr h, c_hssp)
				{
					s = h->seq;
					ba::erase_all(s, "-");
					
					AddSeq(&rs, const_cast<char*>(h->id.c_str()), const_cast<char*>(s.c_str()));
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
				
				c_hssp.clear();
				for (int i = 1; i < msa->nseqs; ++i)
				{
					hit_ptr hit = CreateHit(msa->sqinfo[i].name, msa->seq[0], msa->seq[i]);
					
//					assert(hit->saln.length() == seq.length());
					if (hit->saln.length() < seq.length())
						hit->saln.append(seq.length() - hit->saln.length(), ' ');
					
					if (hit->IdentityAboveThreshold())
					{
						// clean up msa
						uint32 j = 1;
						char* seq = msa->seq[i];
						for (char* si = seq; *si; ++si)
						{
							if (*si == '-')
								continue;
							
							if (j < hit->jfir or j > hit->jlas)
								*si = ' ';
							++j;
						}
						
						c_hssp.push_back(hit);
					}
					else if (VERBOSE)
						cerr << "dropping " << hit->id << endl;
				}
			}
			
			hssp.insert(hssp.end(), c_hssp.begin(), c_hssp.end());

//			cout << endl;
//			for (int i = 0; i < msa->nseqs; ++i)
//			{
//				string id = msa->sqinfo[i].name;
//				id.insert(id.end(), 12 - id.length(), ' ');
//				
//				cout << id << ' ' << msa->seq[i] << endl;
//			}

//			FreeMSeq(&msa);

			uint32 alignmentlength = strlen(msa->seq[0]);
			vector<MResidue*>::const_iterator r = residues.begin();
			for (uint32 i = 0; i < alignmentlength; ++i)
			{
				if (msa->seq[0][i] == '-')
					continue;
				
				assert(r != residues.end());
				assert(kResidueInfo[(*r)->GetType()].code == msa->seq[0][i]);
				
				result.push_back(CreateResidueHInfo(*r, msa, i));
				
				++r;
			}
			
			assert(r == residues.end());
		}
	}
	
	sort(hssp.begin(), hssp.end(), compare_hit());
	uint32 nr = 1;
	foreach (hit_ptr h, hssp)
		h->nr = nr++;

	


//		cout << h << endl;
//	}
//	
//	cout << endl;
//	foreach (hit& h, hssp)
//		cout << h.seq << endl;

	// finally create a HSSP file
//	string hssp;

	using namespace boost::gregorian;
	date today = day_clock::local_day();
	
	// print the header
	os << "HSSP       HOMOLOGY DERIVED SECONDARY STRUCTURE OF PROTEINS , VERSION 2.0d1 2011" << endl
	   << "PDBID      " << inProtein.GetID() << endl
	   << "DATE       file generated on " << to_iso_extended_string(today) << endl
	   << "SEQBASE    " << inDatabank->GetVersion() << endl
	   << "THRESHOLD  according to: t(L)=(290.15 * L ** -0.562) + 5" << endl
	   << "CONTACT    New version by Maarten L. Hekkelman <m.hekkelman@cmbi.ru.nl>" << endl
	   << "HEADER     " + inProtein.GetHeader().substr(10, 40) << endl
	   << "COMPND     " + inProtein.GetCompound().substr(10) << endl
	   << "SOURCE     " + inProtein.GetSource().substr(10) << endl
	   << "AUTHOR     " + inProtein.GetAuthor().substr(10) << endl
//	   << boost::format("SEQLENGTH  %4.4d") % seqlength << endl
//	   << boost::format("NCHAIN     %4.4d chain(s) in %s data set") % nchain % inProtein.GetID() << endl
	   << boost::format("NALIGN     %4.4d") % hssp.size() << endl
	   << endl
	   << "## PROTEINS : EMBL/SWISSPROT identifier and alignment statistics" << endl
	   << "  NR.    ID         STRID   %IDE %WSIM IFIR ILAS JFIR JLAS LALI NGAP LGAP LSEQ2 ACCNUM     PROTEIN" << endl;
	   
	// print the first list
	nr = 1;
	boost::format fmt1("%5.5d : %12.12s%4.4s    %4.2f  %4.2f %4.4d %4.4d %4.4d %4.4d %4.4d %4.4d %4.4d %4.4d  %10.10s %s");
	foreach (hit_ptr h, hssp)
	{
		string id = h->id;
		if (id.length() > 12)
			id.erase(12, string::npos);
		else if (id.length() < 12)
			id.append(12 - id.length(), ' ');
		
		os << fmt1 % nr
				   % id % h->pdb
				   % h->ide % h->wsim % h->ifir % h->ilas % h->jfir % h->jlas % h->lali
				   % h->ngap % h->lgap % h->lseq2
				   % "" % inDatabank->GetMetaData(h->id, "title")
		   << endl;
		
		++nr;
	}

	// print the alignments

	for (uint32 i = 0; i < hssp.size(); i += 70)
	{
		uint32 n = i + 70;
		if (n > hssp.size())
			n = hssp.size();
		
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

		foreach (res_ptr ri, result)
		{
			string aln;
			for (uint32 j = i; j < n; ++j)
				aln += ri->msa->seq[j][ri->pos];
			
			string dssp = ResidueToDSSPLine(inProtein, *ri->res).substr(0, 39);
			
			uint32 nocc = 1, var = 0;
			
			os << ' ' << dssp << boost::format("%4.4d %4.4d  ") % nocc % var << aln << endl;
		}

//		uint32 seqNo = 1;
//		
//		foreach (const MChain* chain, inProtein.GetChains())
//		{
//			foreach (const MResidue* residue, chain->GetResidues())
//			{
//				string aln;
//				for (uint32 j = i; j < n; ++j)
//					aln += hssp[j]->saln[seqNo - 1];
//				
//				const MResidue& r = *residue;
//
//				string dssp = ResidueToDSSPLine(inProtein, *residue).substr(0, 40);
//				
//				uint32 nocc = 1, var = 0;
//				
//				os << ' ' << dssp << boost::format("%4.4d %4.4d  ") % nocc % var << aln << endl;
//				
//				++seqNo;
//			}
//		}		

	}
	
	// ## SEQUENCE PROFILE AND ENTROPY
	os << "## SEQUENCE PROFILE AND ENTROPY" << endl
	   << " SeqNo PDBNo   V   L   I   M   F   W   Y   G   A   P   S   T   C   H   R   K   Q   E   N   D  NOCC NDEL NINS ENTROPY RELENT WEIGHT" << endl;
	
	foreach (res_ptr r, result)
	{
		os << ResidueToDSSPLine(inProtein, *r->res).substr(0, 12);
		
		for (uint32 i = 0; i < 20; ++i)
		{
			os << boost::format("%4.4d") % r->dist[i];
		}
		
		os << "  " << boost::format("%4.4d %4.4d %4.4d") % r->nocc % r->ndel % r->nins << endl;
	}
	
	os << "//" << endl;
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
