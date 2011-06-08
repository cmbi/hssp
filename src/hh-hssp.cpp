//  Copyright Maarten L. Hekkelman, Radboud University 2011.
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
int MAX_HITS = 250;

//namespace hh2
//{
//	
//struct Hit
//{
//	string		seq;
//	uint32		ifir, ilas, jfir, jlas;
//	
//	
//};
//
//class MSA
//{
//  public:
//				MSA();
//				~MSA();
//
//	void		Build(const string& id, const string& seq, CDatabankPtr db);
//
//  private:
//	
//	void		AppendSequence(const string& id, const string& seq);
//
//	mseq_t*		mSA;
//	opts_t		mOpts;
//};
//
//MSA::MSA()
//	: mSA(nil)
//{
//	SetDefaultAlnOpts(&mOpts);
//}
//
//MSA::~MSA()
//{
//	FreeMSeq(&mSA);
//}
//
//void MSA::Build(const string& id, const string& seq, CDatabankPtr inDatabank)
//{
//	NewMSeq(&mSA);
////	AppendSequence(id, seq);
//	AddSeq(&mSA, const_cast<char*>(id.c_str()), const_cast<char*>(seq.c_str()));
//
//	// blast parameters
//	float expect = 1.0;
//	bool filter = true, gapped = true;
//	int wordsize = 3, gapOpen = 11, gapExtend = 1, maxhits = MAX_HITS;
//	string matrix = "BLOSUM62";
//
//	CDbAllDocIterator data(inDatabank.get());
//	CBlast blast(seq, matrix, wordsize, expect, filter, gapped, gapOpen, gapExtend, maxhits);
//	
//	if (blast.Find(*inDatabank, data, nrOfThreads))
//	{
//		CBlastHitList hits(blast.Hits());
//		
//		foreach (const CBlastHit& hit, hits)
//		{
//			const CBlastHsp hsp(hit.Hsps().front());
//			
//			uint32 docNr = hit.DocumentNr();
//			
//			string seq;
//			inDatabank->GetSequence(docNr, inDatabank->GetSequenceNr(docNr, hit.SequenceID()), seq);
//			AppendSequence(hit.DocumentID(), seq);
//		}
//	}
//}
//
//void MSA::AppendSequence(const string& id, const string& seq)
//{
//	mseq_t* nmsa;
//	NewMSeq(&nmsa);
//	AddSeq(&nmsa, const_cast<char*>(id.c_str()), const_cast<char*>(seq.c_str()));
//	if (Align(nmsa, mSA, &mOpts))
//		THROW(("Failure creating multiple sequence alignment"));
//	
//	assert(nmsa->nseqs == mSA->nseqs + 1);
//	FreeMSeq(&mSA);
//	mSA = nmsa;
//}
//	
//}



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

	bool		operator<(const hit& rhs) const 	{ return ide > rhs.ide or (ide == rhs.ide and lali > rhs.lali); }
	
	bool		IdentityAboveThreshold() const;
};

typedef shared_ptr<hit> hit_ptr;

struct ResidueHInfo
{
	char			letter;
	char			chain;
	string			dssp;
	uint32			seqNr;
	uint32			pos;
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

ostream& operator<<(ostream& os, const hit& h)
{
	static boost::format fmt("%5.5d : %12.12s%4.4s    %4.2f  %4.2f %4.4d %4.4d %4.4d %4.4d %4.4d %4.4d %4.4d %4.4d ");
	
	os << fmt % h.nr % h.id % h.pdb
			  % h.ide % h.wsim % h.ifir % h.ilas % h.jfir % h.jlas % h.lali
			  % h.ngap % h.lgap % h.lseq2;
	
	return os;
}

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

	if (VERBOSE and kHomologyThreshold[l] >= ide)
		cerr << "dropping " << id << " because identity " << ide << " is below threshold " << kHomologyThreshold[l] << endl;

	return kHomologyThreshold[l] < ide;
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

	h.lali = ss.length();
	
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
		}
	}

	h.ide = float(h.identical) / float(h.lali);
	h.wsim = float(h.similar) / float(h.lali);
	
	return result;
}

res_ptr CreateResidueHInfo(char a, uint32 nr, vector<hit_ptr>& hits, uint32 pos)
{
	res_ptr r(new ResidueHInfo);
	
	const string kIX("VLIMFWYGAPSTCHRKQEND");
	
	r->letter = a;
	r->chain = 'A';
	r->seqNr = nr;
	r->pos = pos;
	r->nocc = 1;
	r->ndel = r->nins = r->var = 0;
	fill(r->dist, r->dist + 20, 0);
	
	string::size_type ix = kIX.find(a);
	assert(ix != string::npos);
	r->dist[ix] = 1;
	
	foreach (hit_ptr hit, hits)
	{
		if (hit->ifir > nr or hit->ilas < nr)
			continue;
		
		ix = kIX.find(hit->seq[pos]);
		if (ix != string::npos)
		{
			++r->nocc;
			r->dist[ix] += 1;
		}
	}
	
	for (uint32 a = 0; a < 20; ++a)
		r->dist[a] = uint32(((100.0 * r->dist[a]) / r->nocc) + 0.5);

	return r;
}

struct MSAInfo
{
	string				seq;
	vector<char>		chainNames;
	set<hit_ptr>		hits;
	set<res_ptr>		residues;

						MSAInfo(const string& seq, char chainName,
							vector<hit_ptr>& h, vector<res_ptr>& r)
							: seq(seq)
						{
							assert(not h.empty());
							assert(not r.empty());
							
							chainNames.push_back(chainName);
							hits.insert(h.begin(), h.end());
							residues.insert(r.begin(), r.end());
						}
};

char SelectAlignedLetter(const vector<MSAInfo>& msas, hit_ptr hit, res_ptr res)
{
	char result = ' ';
	
	foreach (const MSAInfo& msa, msas)
	{
		if (msa.hits.count(hit) and msa.residues.count(res) and
			hit->ifir <= res->seqNr and hit->ilas >= res->seqNr)
		{
			result = hit->seq[res->pos];
			break;
		}
	}
	
	return result;
}

void ChainToHits(
	CDatabankPtr		inDatabank,
	const string&		seq,
	const string&		seqId,
	opts_t&				coo,
	vector<hit_ptr>&	hssp,
	vector<res_ptr>&	residues)
{
	// blast parameters
	float expect = 1.0;
	bool filter = true, gapped = true;
	int wordsize = 3, gapOpen = 11, gapExtend = 1, maxhits = MAX_HITS;
	string matrix = "BLOSUM62";

	CDbAllDocIterator data(inDatabank.get());
	CBlast blast(seq, matrix, wordsize, expect, filter, gapped, gapOpen, gapExtend, maxhits);
	
	if (blast.Find(*inDatabank, data, nrOfThreads))
	{
		CBlastHitList hits(blast.Hits());
		
		// now create the alignment using clustalo
		mseq_t* msa;
		NewMSeq(&msa);

		AddSeq(&msa, const_cast<char*>(seqId.c_str()), const_cast<char*>(seq.c_str()));
		
		foreach (const CBlastHit& hit, hits)
		{
			const CBlastHsp hsp(hit.Hsps().front());
			
			string s, id;

			inDatabank->GetSequence(hit.DocumentNr(),
				inDatabank->GetSequenceNr(hit.DocumentNr(), hit.SequenceID()),
				s);
			id = hit.DocumentID();
			
			AddSeq(&msa, const_cast<char*>(id.c_str()), const_cast<char*>(s.c_str()));
		}
		
		msa->seqtype = SEQTYPE_PROTEIN;
		msa->aligned = false;
		
		if (Align(msa, nil, &coo))
		{
			FreeMSeq(&msa);
			throw mas_exception("Fatal error creating alignment");
		}
		
		for (int i = 1; i < msa->nseqs; ++i)
		{
			if (msa->sqinfo[i].name == seqId)
				continue;
			
			hit_ptr hit = CreateHit(msa->sqinfo[i].name, msa->seq[0], msa->seq[i]);

			if (hit->IdentityAboveThreshold())
				hssp.push_back(hit);
		}
		
		if (hssp.size() + 1 < msa->nseqs)	// repeat alignment with the new, smaller set of remaining hits
		{
			mseq_t* rs;
			NewMSeq(&rs);
			
			string s = seq;
			ba::erase_all(s, "-");
			
			AddSeq(&rs, const_cast<char*>(seqId.c_str()), const_cast<char*>(s.c_str()));

			foreach (hit_ptr h, hssp)
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
			
			hssp.clear();
			for (int i = 1; i < msa->nseqs; ++i)
			{
				hit_ptr hit = CreateHit(msa->sqinfo[i].name, msa->seq[0], msa->seq[i]);
				
				if (hit->IdentityAboveThreshold())
					hssp.push_back(hit);
			}
		}

		uint32 seqNr = 1;
		for (char* si = msa->seq[0]; *si; ++si)
		{
			if (*si != '-' and *si != ' ')
			{
				residues.push_back(CreateResidueHInfo(*si, seqNr, hssp, si - msa->seq[0]));
				++seqNr;
			}
		}
		
		FreeMSeq(&msa);
	}
}

void CreateHSSP(CDatabankPtr inDatabank, MProtein& inProtein, opts_t& coo, ostream& os)
{
	uint32 nchain = 0, kchain = 0, seqlength = 0;

	vector<hit_ptr> hits;
	vector<res_ptr> res;
	vector<MSAInfo> msas;
	
	foreach (const MChain* chain, inProtein.GetChains())
	{
		const vector<MResidue*>& residues(chain->GetResidues());
		
		string seq;
		chain->GetSequence(seq);
		
		++nchain;
		
		bool newSequence = true;
		foreach (MSAInfo& msa, msas)
		{
			if (msa.seq == seq)
			{
				msa.chainNames.push_back(chain->GetChainID());
				newSequence = false;
			}
		}

		if (not newSequence)
			continue;
		
		++kchain;
		seqlength += seq.length();

		vector<hit_ptr> c_hits;
		
		ChainToHits(inDatabank, seq, inProtein.GetID(), coo, c_hits, res);
		if (not c_hits.empty() and not res.empty())
		{
			assert(res.size() == residues.size());
			for (uint32 i = 0; i < res.size(); ++i)
			{
				assert(kResidueInfo[residues[i]->GetType()].code == res[i]->letter);
				assert(residues[i]->GetSeqNumber() == res[i]->seqNr);
				res[i]->chain = chain->GetChainID();
				res[i]->dssp = ResidueToDSSPLine(inProtein, *residues[i]).substr(0, 39);
			}
			
			msas.push_back(MSAInfo(seq, chain->GetChainID(), c_hits, res));
			hits.insert(hits.end(), c_hits.begin(), c_hits.end());
		}
	}
	
	sort(hits.begin(), hits.end(), compare_hit());
	uint32 nr = 1;
	foreach (hit_ptr h, hits)
		h->nr = nr++;

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
	   << boost::format("SEQLENGTH  %4.4d") % seqlength << endl
	   << boost::format("NCHAIN     %4.4d chain(s) in %s data set") % nchain % inProtein.GetID() << endl;
	
	if (kchain != nchain)
	{
		os << boost::format("KCHAIN     %4.4d chain(s) used here ; chains(s) : ") % kchain << endl;
	}
	
	os << boost::format("NALIGN     %4.4d") % hits.size() << endl
	   << endl
	   << "## PROTEINS : EMBL/SWISSPROT identifier and alignment statistics" << endl
	   << "  NR.    ID         STRID   %IDE %WSIM IFIR ILAS JFIR JLAS LALI NGAP LGAP LSEQ2 ACCNUM     PROTEIN" << endl;
	   
	// print the first list
	nr = 1;
	boost::format fmt1("%5.5d : %12.12s%4.4s    %4.2f  %4.2f %4.4d %4.4d %4.4d %4.4d %4.4d %4.4d %4.4d %4.4d  %10.10s %s");
	foreach (hit_ptr h, hits)
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
	for (uint32 i = 0; i < hits.size(); i += 70)
	{
		uint32 n = i + 70;
		if (n > hits.size())
			n = hits.size();
		
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

		foreach (res_ptr ri, res)
		{
			string aln;
			
			for (uint32 j = i; j < n; ++j)
				aln += SelectAlignedLetter(msas, hits[j], ri);
			
			os << ' ' << ri->dssp << boost::format("%4.4d %4.4d  ") % ri->nocc % ri->var << aln << endl;
		}
	}
	
	// ## SEQUENCE PROFILE AND ENTROPY
	os << "## SEQUENCE PROFILE AND ENTROPY" << endl
	   << " SeqNo PDBNo   V   L   I   M   F   W   Y   G   A   P   S   T   C   H   R   K   Q   E   N   D  NOCC NDEL NINS ENTROPY RELENT WEIGHT" << endl;
	
	uint32 seqNr = 1;
	foreach (res_ptr r, res)
	{
		os << boost::format(" %4.4d %4.4d %c") % seqNr % seqNr % r->chain;
		++seqNr;

		for (uint32 i = 0; i < 20; ++i)
			os << boost::format("%4.4d") % r->dist[i];
		
		os << "  " << boost::format("%4.4d %4.4d %4.4d") % r->nocc % r->ndel % r->nins << endl;
	}
	
	os << "//" << endl;
}

void CreateHSSP(CDatabankPtr inDatabank, const string& inProtein, opts_t& coo, ostream& os)
{
	vector<hit_ptr> hits;
	vector<res_ptr> res;

	string seq = inProtein;
	ba::erase_all(seq, "\r\n");
	ba::erase_all(seq, "\n");

	ChainToHits(inDatabank, seq, "UNKN", coo, hits, res);

	vector<MSAInfo> msas;
	msas.push_back(MSAInfo(seq, 'A', hits, res));
	
	sort(hits.begin(), hits.end(), compare_hit());
	uint32 nr = 1;
	foreach (hit_ptr h, hits)
		h->nr = nr++;
	
	using namespace boost::gregorian;
	date today = day_clock::local_day();
	
	// print the header
	os << "HSSP       HOMOLOGY DERIVED SECONDARY STRUCTURE OF PROTEINS , VERSION 2.0d1 2011" << endl
	   << "DATE       file generated on " << to_iso_extended_string(today) << endl
	   << "SEQBASE    " << inDatabank->GetVersion() << endl
	   << "THRESHOLD  according to: t(L)=(290.15 * L ** -0.562) + 5" << endl
	   << "CONTACT    New version by Maarten L. Hekkelman <m.hekkelman@cmbi.ru.nl>" << endl
	   << boost::format("SEQLENGTH  %4.4d") % seq.length() << endl
	   << "NCHAIN     1 chain(s) in data set" << endl
	   << boost::format("NALIGN     %4.4d") % hits.size() << endl
	   << endl
	   << "## PROTEINS : EMBL/SWISSPROT identifier and alignment statistics" << endl
	   << "  NR.    ID         STRID   %IDE %WSIM IFIR ILAS JFIR JLAS LALI NGAP LGAP LSEQ2 ACCNUM     PROTEIN" << endl;
	   
	// print the first list
	nr = 1;
	boost::format fmt1("%5.5d : %12.12s%4.4s    %4.2f  %4.2f %4.4d %4.4d %4.4d %4.4d %4.4d %4.4d %4.4d %4.4d  %10.10s %s");
	foreach (hit_ptr h, hits)
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
	for (uint32 i = 0; i < hits.size(); i += 70)
	{
		uint32 n = i + 70;
		if (n > hits.size())
			n = hits.size();
		
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

		uint32 seqNr = 1;
		foreach (res_ptr ri, res)
		{
			string aln;
			
			for (uint32 j = i; j < n; ++j)
				aln += SelectAlignedLetter(msas, hits[j], ri);
			
			os << boost::format("  %4.4d %4.4d A %c                         %4.4d %4.4d  ")
				% seqNr % seqNr % ri->letter % ri->nocc % ri->var << aln << endl;

			++seqNr;
		}
	}
	
	// ## SEQUENCE PROFILE AND ENTROPY
	os << "## SEQUENCE PROFILE AND ENTROPY" << endl
	   << " SeqNo PDBNo   V   L   I   M   F   W   Y   G   A   P   S   T   C   H   R   K   Q   E   N   D  NOCC NDEL NINS ENTROPY RELENT WEIGHT" << endl;
	
	uint32 seqNr = 1;
	foreach (res_ptr r, res)
	{
		os << boost::format(" %4.4d %4.4d A") % seqNr % seqNr;
		++seqNr;
		
		for (uint32 i = 0; i < 20; ++i)
			os << boost::format("%4.4d") % r->dist[i];
		
		os << "  " << boost::format("%4.4d %4.4d %4.4d") % r->nocc % r->ndel % r->nins << endl;
	}
	
	os << "//" << endl;
}

namespace hh
{

void Init()
{
	static bool sInited = false;
	if (not sInited)
	{
		sInited = true;

	    LogDefaultSetup(&rLog);

#if 1 // DEBUG
	    rLog.iLogLevelEnabled = LOG_VERBOSE;
	    VERBOSE = 1;
#endif

		InitClustalOmega(BLAST_THREADS);
	}
}

void CreateHSSP(
	CDatabankPtr				inDatabank,
	MProtein&					inProtein,
	std::ostream&				outHSSP)
{
	Init();

	opts_t coo;
	SetDefaultAlnOpts(&coo);
	
	CreateHSSP(inDatabank, inProtein, coo, outHSSP);
}

void CreateHSSP(
	CDatabankPtr				inDatabank,
	const std::string&			inProtein,
	std::ostream&				outHSSP)
{
	Init();
	
	opts_t coo;
	SetDefaultAlnOpts(&coo);
	
	CreateHSSP(inDatabank, inProtein, coo, outHSSP);
}

}
