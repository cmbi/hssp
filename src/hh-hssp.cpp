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
#include <boost/iostreams/device/file_descriptor.hpp>
#include <boost/foreach.hpp>
#define foreach BOOST_FOREACH
#include <boost/program_options.hpp>
#include <boost/ptr_container/ptr_vector.hpp>
#include <boost/date_time/gregorian/gregorian.hpp>
#include <boost/date_time/date_clock_device.hpp>
#include <boost/regex.hpp>

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
#include "hh-hssp.h"

using namespace std;
namespace ba = boost::algorithm;
namespace io = boost::iostreams;
namespace po = boost::program_options;

int	nrOfThreads = boost::thread::hardware_concurrency(),
	MAX_HITS = 250;

namespace hh
{

void AlignWithClustalOmega(const string& inClustalO, mseq& seqs)
{
	HUuid uuid;
	
	fs::path rundir("/tmp/hssp-2/");
	rundir /= boost::lexical_cast<string>(uuid);
	fs::create_directories(rundir);
	
	// start a clustalo
	int ifd[2];
	int ofd[2];
	
	if (pipe(ifd) != 0)
		THROW(("Failed to create pipe: %s", strerror(errno)));

	if (pipe(ofd) != 0)
	{
		close(ifd[0]);
		close(ifd[1]);
		THROW(("Failed to create pipe: %s", strerror(errno)));
	}
	
	int pid = fork();
	
	if (pid == -1)
	{
		close(ifd[0]);
		close(ifd[1]);
		close(ofd[0]);
		close(ofd[1]);
		
		THROW(("fork failed: %s", strerror(errno)));
	}
	
	if (pid == 0)	// the child process (will be clustalo)
	{
		fs::current_path(rundir);
		
		setpgid(0, 0);
		
		dup2(ifd[0], STDIN_FILENO);
		close(ifd[0]);
		close(ifd[1]);
		
		int fd = open("clustalo.log", O_CREAT | O_RDWR | O_APPEND, 0666);
//		dup2(fd, STDOUT_FILENO);
		dup2(fd, STDERR_FILENO);
		
		dup2(ofd[1], STDOUT_FILENO);
		close(ofd[0]);
		close(ofd[1]);
		
		char* argv[] = {
			const_cast<char*>(inClustalO.c_str()),
			const_cast<char*>("-i"),
			const_cast<char*>("-"),
			nil
		};
		
		(void)execve(inClustalO.c_str(), argv, environ);
	}

	close(ofd[1]);
	close(ifd[0]);
	int fd = ifd[1];
	
	// ClustalO is running and waiting for a FastA file on stdin
	foreach (const seq& s, seqs)
	{
		stringstream ss;
		ss << '>' << s.m_id << endl;
		
		for (uint32 o = 0; o < s.m_seq.length(); o += 80)
		{
			uint32 n = s.m_seq.length() - o;
			if (n > 80)
				n = 80;
			
			ss << s.m_seq.substr(o, n) << endl;
		}

		WriteToFD(fd, ss.str());
	}

	close(fd);	// Closing fd will tell ClustalO we're done feeding input
	
	// now read in the result
	
	mseq result;
	
	io::filtering_istream is;
	is.push(io::file_descriptor_source(ofd[0]));
	
	for (;;)
	{
		string line;
		getline(is, line);

		if (line.empty() and is.eof())
			break;

		if (ba::starts_with(line, ">>"))
		{
			cerr << "Invalid clustalo output, crash?" << endl;
			break;
		}
		
		if (ba::starts_with(line, ">"))
		{
			result.push_back(seq());
			result.back().m_id = line.substr(1);
			continue;
		}
		
		if (result.empty())
		{
			cerr << "Invalid clustalo output" << endl;
			break;
		}
		
		result.back().m_seq.append(line);
	}
	
	close(ofd[0]);

	int status;
	waitpid(pid, &status, 0);
	
	if (status != 0)
		THROW(("clustalo exited with status %d", status));
	
	swap(seqs, result);
}
	
struct insertion
{
	uint32		ipos, jpos;
	string		seq;
};

struct hit
{
				hit(const string& id, const string& seq, const string& desc) : nr(0), id(id), seq(seq), desc(desc), alwaysSelect(false) {}

	uint32		nr;
	string		id, acc, seq, desc, pdb;
	uint32		ifir, ilas, jfir, jlas, lali, ngap, lgap, lseq2;
	float		ide, wsim;
	uint32		identical, similar;
	vector<insertion>
				insertions;
	bool		alwaysSelect;

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

hit_ptr CreateHit(const string& id, const string& q, const string& s, const string& d)
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

	result.reset(new hit(id, s, d));
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
		if (not hit->alwaysSelect and (hit->ifir > nr or hit->ilas < nr))
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
			(hit->alwaysSelect or (hit->ifir <= res->seqNr and hit->ilas >= res->seqNr)))
		{
			result = hit->seq[res->pos];
			break;
		}
	}
	
	return result;
}

void ChainToHits(
	CDatabankPtr		inDatabank,
	const string&		inClustalO,
	const string&		seq,
	const string&		seqId,
	vector<hit_ptr>&	hssp,
	vector<res_ptr>&	residues)
{
	// blast parameters
	float expect = 1.0;
	bool filter = true, gapped = true;
	int wordsize = 3, gapOpen = 11, gapExtend = 1, maxhits = MAX_HITS;
	string matrix = "BLOSUM62";

	CBlastResult* result = PerformBlastSearch(*inDatabank, seq, matrix, wordsize, expect, filter, gapped, gapOpen, gapExtend, maxhits);
	
	if (result != nil)
	{
		CBlastHitList& hits = result->hits;
		
		// now create the alignment using clustalo
		hh::mseq msa;
		msa.push_back(hh::seq(seqId, seq));
		
		foreach (const CBlastHit& hit, hits)
		{
			const CBlastHsp hsp(hit.hsps.front());
			
			string s, id;

			inDatabank->GetSequence(hit.documentNr,
				inDatabank->GetSequenceNr(hit.documentNr, hit.sequenceID),
				s);
			id = hit.documentID;
			
			msa.push_back(hh::seq(id, s));
		}
		
		hh::AlignWithClustalOmega(inClustalO, msa);
		
		for (int i = 1; i < msa.size(); ++i)
		{
			hit_ptr hit = CreateHit(msa[i].m_id, msa[0].m_seq, msa[i].m_seq, msa[i].m_desc);

			if (hit->IdentityAboveThreshold())
				hssp.push_back(hit);
		}
		
		if (hssp.size() + 1 < msa.size())	// repeat alignment with the new, smaller set of remaining hits
		{
			msa.clear();

			msa.push_back(hh::seq(seqId, seq));

			foreach (hit_ptr h, hssp)
			{
				string s = h->seq;
				ba::erase_all(s, "-");
				
				msa.push_back(hh::seq(h->id, s));
			}
			
			hh::AlignWithClustalOmega(inClustalO, msa);
			hssp.clear();
			
			for (int i = 1; i < msa.size(); ++i)
			{
				hit_ptr hit = CreateHit(msa[i].m_id, msa[0].m_seq, msa[i].m_seq, msa[i].m_desc);
	
				if (hit->IdentityAboveThreshold())
					hssp.push_back(hit);
			}
		}

		uint32 seqNr = 1;
		string s = msa.front().m_seq;
		for (uint32 i = 0; i < s.length(); ++i)
		{
			if (s[i] != '-' and s[i] != ' ')
			{
				residues.push_back(CreateResidueHInfo(s[i], seqNr, hssp, i));
				++seqNr;
			}
		}
	}
}

void CreateHSSPOutput(
	const MProtein&		inProtein,
	const string&		inDatabankVersion,
	uint32				inSeqLength,
	uint32				inKChain,
	uint32				inNChain,
	vector<hit_ptr>&	hits,
	vector<res_ptr>&	res,
	vector<MSAInfo>&	msas,
	ostream&			os)
{
	using namespace boost::gregorian;
	date today = day_clock::local_day();
	
	// print the header
	os << "HSSP       HOMOLOGY DERIVED SECONDARY STRUCTURE OF PROTEINS , VERSION 2.0d1 2011" << endl
	   << "PDBID      " << inProtein.GetID() << endl
	   << "DATE       file generated on " << to_iso_extended_string(today) << endl
	   << "SEQBASE    " << inDatabankVersion << endl
	   << "THRESHOLD  according to: t(L)=(290.15 * L ** -0.562) + 5" << endl
	   << "CONTACT    New version by Maarten L. Hekkelman <m.hekkelman@cmbi.ru.nl>" << endl
	   << "HEADER     " + inProtein.GetHeader().substr(10, 40) << endl
	   << "COMPND     " + inProtein.GetCompound().substr(10) << endl
	   << "SOURCE     " + inProtein.GetSource().substr(10) << endl
	   << "AUTHOR     " + inProtein.GetAuthor().substr(10) << endl
	   << boost::format("SEQLENGTH  %4.4d") % inSeqLength << endl
	   << boost::format("NCHAIN     %4.4d chain(s) in %s data set") % inNChain % inProtein.GetID() << endl;
	
	if (inKChain != inNChain)
	{
		os << boost::format("KCHAIN     %4.4d chain(s) used here ; chains(s) : ") % inKChain << endl;
	}
	
	os << boost::format("NALIGN     %4.4d") % hits.size() << endl
	   << endl
	   << "## PROTEINS : EMBL/SWISSPROT identifier and alignment statistics" << endl
	   << "  NR.    ID         STRID   %IDE %WSIM IFIR ILAS JFIR JLAS LALI NGAP LGAP LSEQ2 ACCNUM     PROTEIN" << endl;
	   
	// print the first list
	uint32 nr = 1;
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
				   % "" % h->desc
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

void CreateHSSP(CDatabankPtr inDatabank, const string& inClustalO, MProtein& inProtein, ostream& os)
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
		
		ChainToHits(inDatabank, inClustalO, seq, inProtein.GetID(), c_hits, res);
		if (not c_hits.empty() and not res.empty())
		{
			assert(res.size() == residues.size());
			for (uint32 i = 0; i < res.size(); ++i)
			{
				assert(kResidueInfo[residues[i]->GetType()].code == res[i]->letter);
				assert(residues[i]->GetSeqNumber() == res[i]->seqNr);
				res[i]->chain = chain->GetChainID();
				res[i]->dssp = ResidueToDSSPLine(*residues[i]).substr(0, 39);
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

void CreateHSSP(CDatabankPtr inDatabank, const string& inClustalO, const string& inProtein, ostream& os)
{
	vector<hit_ptr> hits;
	vector<res_ptr> res;

	string seq = inProtein;
	ba::erase_all(seq, "\r\n");
	ba::erase_all(seq, "\n");

	ChainToHits(inDatabank, inClustalO, seq, "UNKN", hits, res);

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

void CreateHSSPForAlignment(
	const mseq&					inAlignment,
	MProtein&					inProtein,
	char						inChain,
	ostream&					outHSSP)
{
	boost::regex re("([a-zA-Z0-9_]+)/(\\d+)-(\\d+)");

	vector<hit_ptr> hits;
	for (int i = 1; i < inAlignment.size(); ++i)
	{
		hit_ptr hit = CreateHit(inAlignment[i].m_id, inAlignment[0].m_seq, inAlignment[i].m_seq, inAlignment[i].m_desc);

		// parse out the position
		boost::smatch sm;
		if (not boost::regex_match(inAlignment[i].m_id, sm, re))
			throw mas_exception("Alignment ID should contain position");
		
		hit->id = sm.str(1);
		hit->ifir = boost::lexical_cast<uint32>(sm.str(2));
		hit->ilas = boost::lexical_cast<uint32>(sm.str(3));
		hit->alwaysSelect = true;

		for (string::iterator r = hit->seq.begin(); r != hit->seq.end() and *r == '-'; ++r)
			*r = ' ';

		for (string::reverse_iterator r = hit->seq.rbegin(); r != hit->seq.rend() and *r == '-'; ++r)
			*r = ' ';

		if (hit->IdentityAboveThreshold())
			hits.push_back(hit);
	}

	sort(hits.begin(), hits.end(), compare_hit());
	if (hits.size() > 9999)
		hits.erase(hits.begin() + 9999, hits.end());
	
	uint32 nr = 1;
	foreach (hit_ptr h, hits)
		h->nr = nr++;
	
	string seq;
	inProtein.GetChain(inChain).GetSequence(seq);

	const vector<MResidue*>& residues(inProtein.GetChain(inChain).GetResidues());
	vector<res_ptr> res;
	uint32 seqNr = 0;
	string s = inAlignment.front().m_seq;
	for (uint32 i = 0; i < s.length(); ++i)
	{
		if (s[i] != '-' and s[i] != ' ' and s[i] != '.')
		{
			res.push_back(CreateResidueHInfo(s[i], seqNr + 1, hits, i));

			assert(kResidueInfo[residues[seqNr]->GetType()].code == s[i]);
			assert(s[i] == seq[seqNr]);
			res.back()->seqNr = residues[seqNr]->GetSeqNumber();
			res.back()->chain = inChain;
			res.back()->dssp = ResidueToDSSPLine(*residues[seqNr]).substr(0, 39);

			++seqNr;
		}
	}

	assert(res.size() == seq.length());

	vector<MSAInfo> msas;
	msas.push_back(MSAInfo(seq, inChain, hits, res));
	
	CreateHSSPOutput(inProtein, "", seq.length(), 1, 1, hits, res, msas, outHSSP);
}

}
