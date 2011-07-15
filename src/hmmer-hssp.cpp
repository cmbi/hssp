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
#include "hmmer-hssp.h"

using namespace std;
namespace ba = boost::algorithm;
namespace io = boost::iostreams;
namespace po = boost::program_options;

namespace hmmer
{
	
inline bool is_gap(char aa)
{
	return aa == '-' or aa == '~' or aa == '.' or aa == '_';
}

struct seq
{
	std::string		m_id;
	std::string		m_seq;
	
					seq() {}

					seq(const std::string& id, const std::string& seq)
						: m_id(id)
						, m_seq(seq)
					{
					}
};

typedef std::vector<seq> mseq;

// ReadStockholm is a function that reads a multiple sequence alignment from
// a Stockholm formatted file. Restriction is that this Stockholm file has
// a #=GF field at the second line containing the ID of the query used in
// jackhmmer.	

void ReadStockholm(istream& is, mseq& msa)
{
	string line;
	getline(is, line);
	if (line != "# STOCKHOLM 1.0")
		throw mas_exception("Not a stockholm file");

	getline(is, line);
	if (not ba::starts_with(line, "#=GF ID "))
		throw mas_exception("Not a valid stockholm file, missing #=GF ID line");
	
	string id = line.substr(8);

	boost::regex re("(.+?)-i(?:\\d+)$");
	boost::smatch sm;
	if (boost::regex_match(id, sm, re))
		id = sm.str(1);

	map<string,uint32> ix;
	
	// fill the index
	msa.push_back(seq(id, ""));
	ix[id] = 0;
	
	for (;;)
	{
		getline(is, line);
		
		if (line.empty())
		{
			if (is.eof())
				break;
			continue;
		}
		
		if (line == "//")
			break;
		
		if (ba::starts_with(line, "#=GS "))
		{
			uint32 n = msa.size();
			msa.push_back(hmmer::seq());
			
			string id = line.substr(5);
			string::size_type s = id.find("DE ");
			if (s != string::npos)
				id = id.substr(0, s);
			
			ba::trim(id);
			msa[n].m_id = id;
			ix[id] = n;
			continue;
		}
		
		if (line[0] != '#')
		{
			string::size_type s = line.find(' ');
			if (s == string::npos)
				throw mas_exception("Invalid stockholm file");
			
			string id = line.substr(0, s);
			
			while (s < line.length() and line[s] == ' ')
				++s;
			
			string seq = line.substr(s);
			
			map<string,uint32>::iterator i = ix.find(id);
			if (i == ix.end())
			{
				ix.insert(make_pair(id, msa.size()));
				msa.push_back(hmmer::seq());
				msa.back().m_id = id;
				i = ix.find(id);
			}
			
			msa[i->second].m_seq += seq;
		}
	}
}

void RunJackHmmer(const string& seq, uint32 iterations, const fs::path& fastadir, const fs::path& jackhmmer,
	const string& db, mseq& msa)
{
	HUuid uuid;
	
	fs::path rundir("/tmp/hssp-2/");
	rundir /= boost::lexical_cast<string>(uuid);
	fs::create_directories(rundir);
	
	// write fasta file
	
	fs::ofstream input(rundir / "input.fa");
	if (not input.is_open())
		throw mas_exception("Failed to create jackhmmer input file");
		
	input << '>' << "input" << endl;
	for (uint32 o = 0; o < seq.length(); o += 72)
	{
		uint32 k = seq.length() - o;
		if (k > 72)
			k = 72;
		input << seq.substr(o, k) << endl;
	}
	input.close();
	
	// start a jackhmmer
	int pid = fork();
	
	if (pid == -1)
		THROW(("fork failed: %s", strerror(errno)));
	
	if (pid == 0)	// the child process (will be jackhmmer)
	{
		fs::current_path(rundir);
		
		setpgid(0, 0);
		
		int fd = open("jackhmmer.log", O_CREAT | O_RDWR | O_APPEND, 0666);
		dup2(fd, STDOUT_FILENO);
		dup2(fd, STDERR_FILENO);
		
		arg_vector argv(jackhmmer.string());
		
		argv.push("-N", iterations);
		argv.push("--noali");
//		argv.push("-o", "/dev/null");
		argv.push("-A", "output.sto");
		argv.push("input.fa");
		argv.push((fastadir / (db + ".fa")).string());
		
		if (VERBOSE)
			cerr << argv << endl;
		
		(void)execve(jackhmmer.string().c_str(), argv, environ);
		cerr << "Failed to run " << jackhmmer << endl << " err: " << strerror(errno) << endl;
		exit(-1);
	}

	// wait for jackhmmer to finish
	int status;
	waitpid(pid, &status, 0);
	
	if (status != 0)
	{
		if (fs::exists(rundir / "jackhmmer.log"))
		{
			fs::ifstream log(rundir / "jackhmmer.log");
			if (log.is_open())
				io::copy(log, cerr);
		}
		
		THROW(("jackhmmer exited with status %d", status));
	}

	// read in the result
	if (not fs::exists(rundir / "output.sto"))
		THROW(("Output Stockholm file is missing"));
	
	fs::ifstream is(rundir / "output.sto");
	ReadStockholm(is, msa);
	
	if (not VERBOSE)
		fs::remove_all(rundir);
}
	
struct insertion
{
	uint32		ipos, jpos;
	string		seq;
};

struct hit
{
				hit(const mseq& msa, uint32 ix)
					: msa(msa), ix(ix) {}

	const mseq&	msa;
	uint32		nr, ix;
	string		id, acc, desc, pdb;
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

// Create a hit object based on a jackhmmer alignment pair
// first is the original query sequence, with gaps introduced.
// second is the hit sequence.
// Since this is jackhmmer output, we can safely assume the
// alignment does not contain gaps at the start or end of the query.
// (However, being paranoid, I do check this...)
hit_ptr CreateHit(mseq& msa, uint32 qix, uint32 six)
{
	hit_ptr result;
	
	string& q = msa[qix].m_seq;
	string& s = msa[six].m_seq;
	string& id = msa[six].m_id;
	
	if (q.empty() or s.empty())
		THROW(("Invalid (empty) sequence"));
	
	if (is_gap(q[0]) or is_gap(q[q.length() - 1]))
		THROW(("Leading (or trailing) gaps found in query sequence"));
	
	assert(q.length() == s.length());

	result.reset(new hit(msa, six));
	hit& h = *result;

	// parse out the position
	static const boost::regex re("([a-zA-Z0-9_]+)/(\\d+)-(\\d+)");
	boost::smatch sm;

	if (not boost::regex_match(id, sm, re))
		throw mas_exception("Alignment ID should contain position");
	
	h.id = sm.str(1);
	
	h.ifir = 1;
	h.ilas = q.length();

	// jfir/jlas can be taken over from jackhmmer output
	h.jfir = boost::lexical_cast<uint32>(sm.str(2));
	h.jlas = boost::lexical_cast<uint32>(sm.str(3));

	h.lgap = h.ngap = h.identical = h.similar = 0;
	
	string::iterator qb = q.begin(), qe = q.end(), sb = s.begin(), se = s.end();

	while (qb != qe)
	{
		if (is_gap(*sb))
		{
			*sb = ' ';
			++h.ifir;
		}
		else
			break;
		
		++qb;
		++sb;
	}
	
	while (qe != qb and is_gap(*(se - 1)))
	{
		--h.ilas;
		--qe;
		--se;
		*se = ' ';
	}

	h.lali = s.length();
	
	bool gap = true;
	const substitution_matrix m("BLOSUM62");
	for (string::iterator si = sb, qi = qb; si != se; ++si, ++qi)
	{
		if (is_gap(*si) and is_gap(*qi))
		{
			// a common gap
			--h.lali;
		}
		else if (is_gap(*si))
		{
			if (not gap)
				++h.ngap;
			gap = true;
			++h.lgap;
		}
		else if (is_gap(*qi))
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

res_ptr CreateResidueHInfo(char a, vector<hit_ptr>& hits, uint32 pos)
{
	res_ptr r(new ResidueHInfo);
	
	const string kIX("VLIMFWYGAPSTCHRKQEND");
	
	r->letter = a;
	r->chain = 'A';
	r->seqNr = 0;
	r->pos = pos;
	r->nocc = 1;
	r->ndel = r->nins = r->var = 0;
	fill(r->dist, r->dist + 20, 0);
	
	string::size_type ix = kIX.find(a);
	assert(ix != string::npos);
	r->dist[ix] = 1;
	
	foreach (hit_ptr hit, hits)
	{
		ix = kIX.find(hit->msa[hit->ix].m_seq[pos]);
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
	set<hit_ptr>		hits;
	set<res_ptr>		residues;

						MSAInfo(const string& seq, vector<hit_ptr>& h, vector<res_ptr>& r)
							: seq(seq)
						{
							assert(not h.empty());
							assert(not r.empty());

							hits.insert(h.begin(), h.end());
							residues.insert(r.begin(), r.end());
						}
};

char SelectAlignedLetter(const vector<MSAInfo>& msas, hit_ptr hit, res_ptr res)
{
	char result = ' ';
	
	foreach (const MSAInfo& msa, msas)
	{
		if (msa.hits.count(hit) and msa.residues.count(res))
		{
			uint32 p = res->pos;
			uint32 i = hit->ix;
			
			const string& q = hit->msa[0].m_seq;
			const string& s = hit->msa[i].m_seq;
			
			result = s[p];

			if (not is_gap(result) and
				((p > 0 and is_gap(q[p - 1]) and not is_gap(s[p - 1])) or
				 (p + 1 < q.length() and is_gap(q[p + 1]) and not is_gap(s[p + 1]))))
				result = tolower(result);
			
			break;
		}
	}
	
	return result;
}

void CreateHSSPOutput(
	const MProtein&		inProtein,
	const string&		inDatabankVersion,
	uint32				inSeqLength,
	uint32				inNChain,
	uint32				inKChain,
	const string&		inUsedChains,
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
		os << boost::format("KCHAIN     %4.4d chain(s) used here ; chains(s) : ") % inKChain << inUsedChains << endl;
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
		
		string acc = h->acc;
		if (acc.length() > 10)
			acc.erase(10, string::npos);
		else if (acc.length() < 10)
			acc.append(10 - acc.length(), ' ');
		
		os << fmt1 % nr
				   % id % h->pdb
				   % h->ide % h->wsim % h->ifir % h->ilas % h->jfir % h->jlas % h->lali
				   % h->ngap % h->lgap % h->lseq2
				   % acc % h->desc
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

		res_ptr last;
		uint32 seqNr = 1;
		foreach (res_ptr ri, res)
		{
			if (last and last->seqNr + 1 != ri->seqNr)
			{
				os << boost::format(" %5.5d        !  !           0   0    0    0    0") % seqNr << endl;
				++seqNr;
			}

			last = ri;
			++seqNr;
			
			string aln;
			
			for (uint32 j = i; j < n; ++j)
				aln += SelectAlignedLetter(msas, hits[j], ri);
			
			os << ' ' << boost::format("%5.5d%s%4.4d %4.4d  ") % ri->seqNr % ri->dssp % ri->nocc % ri->var << aln << endl;
		}
	}
	
	// ## SEQUENCE PROFILE AND ENTROPY
	os << "## SEQUENCE PROFILE AND ENTROPY" << endl
	   << " SeqNo PDBNo   V   L   I   M   F   W   Y   G   A   P   S   T   C   H   R   K   Q   E   N   D  NOCC NDEL NINS ENTROPY RELENT WEIGHT" << endl;
	
	uint32 seqNr = 1;
	res_ptr last;
	foreach (res_ptr r, res)
	{
		if (last and last->seqNr + 1 != r->seqNr)
		{
			os << boost::format("%5.5d          0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0     0    0    0   0.000      0  1.00")
				% seqNr << endl;
			++seqNr;
		}
		last = r;

		os << boost::format(" %4.4d %4.4d %c") % seqNr % seqNr % r->chain;
		++seqNr;

		for (uint32 i = 0; i < 20; ++i)
			os << boost::format("%4.4d") % r->dist[i];
		
		os << "  " << boost::format("%4.4d %4.4d %4.4d") % r->nocc % r->ndel % r->nins << endl;
	}
	
	os << "//" << endl;
}

void ChainToHits(CDatabankPtr inDatabank, mseq& msa,
	vector<hit_ptr>& hits, const vector<MResidue*>& residues, vector<res_ptr>& res)
{
	for (uint32 i = 1; i < msa.size(); ++i)
	{
		hit_ptr hit = CreateHit(msa, 0, i);

//		for (string::iterator r = hit->seq.begin(); r != hit->seq.end() and *r == '-'; ++r)
//			*r = ' ';
//
//		for (string::reverse_iterator r = hit->seq.rbegin(); r != hit->seq.rend() and *r == '-'; ++r)
//			*r = ' ';

		if (hit->IdentityAboveThreshold())
		{
			uint32 docNr = inDatabank->GetDocumentNr(hit->id);
			
			hit->desc = inDatabank->GetMetaData(docNr, "title");
			try
			{
				hit->acc = inDatabank->GetMetaData(docNr, "acc");
			}
			catch (...) {}
			hit->lseq2 = inDatabank->GetSequence(docNr, 0).length();
			
			hits.push_back(hit);
		}
	}

	const string& s = msa.front().m_seq;
	for (uint32 i = 0; i < s.length(); ++i)
	{
		if (s[i] != '-' and s[i] != ' ')
			res.push_back(CreateResidueHInfo(s[i], hits, i));
	}
}

// Find the minimal set of overlapping sequences
// Only search fully contained subsequences, no idea what to do with
// sequences that overlap and each have a tail. What residue number to use in that case? What chain ID?
void ClusterSequences(vector<string>& s, vector<uint32>& ix)
{
	for (;;)
	{
		bool found = false;
		for (uint32 i = 0; not found and i < s.size() - 1; ++i)
		{
			for (uint32 j = i + 1; not found and j < s.size(); ++j)
			{
				string& a = s[i];
				string& b = s[j];

				if (a.empty() or b.empty())
					continue;

				if (ba::contains(a, b)) // j fully contained in i
				{
					s[j].clear();
					ix[j] = i;
					found = true;
				}
				else if (ba::contains(b, a)) // i fully contained in j
				{
					s[i].clear();
					ix[i] = j;
					found = true;
				}
			}
		}
		
		if (not found)
			break;
	}
}

void CreateHSSP(
	const std::string&			inProtein,
	std::ostream&				outHSSP)
{
}

void CreateHSSP(
	CDatabankPtr				inDatabank,
	MProtein&					inProtein,
	const fs::path&				inFastaDir,
	const fs::path&				inJackHmmer,
	uint32						inIterations,
	uint32						inMinSeqLength,
	ostream&					outHSSP)
{
	uint32 seqlength = 0;

	vector<hit_ptr> hits;
	vector<res_ptr> res;

	// construct a set of unique sequences, containing only the largest ones in case of overlap
	vector<string> seqset;
	vector<uint32> ix;
	vector<const MChain*> chains;
	
	foreach (const MChain* chain, inProtein.GetChains())
	{
		string seq;
		chain->GetSequence(seq);
		
		if (seq.length() < inMinSeqLength)
			continue;
		
		chains.push_back(chain);
		seqset.push_back(seq);
		ix.push_back(ix.size());
	}
	
	if (seqset.empty())
		THROW(("Not enough sequences in DSSP file of length %d", inMinSeqLength));

	if (seqset.size() > 1)
		ClusterSequences(seqset, ix);

	vector<mseq> alignments;
	for (uint32 i = 0; i < ix.size(); ++i)
	{
		alignments.push_back(mseq());
		if (ix[i] != i)	// if remapped (double) skip
			continue;
		RunJackHmmer(seqset[i], inIterations, inFastaDir, inJackHmmer, inDatabank->GetID(), alignments.back());
	}
	
	vector<MSAInfo> msas;
	uint32 seqNr = 1;
	for (uint32 i = 0; i < chains.size(); ++i)
	{
		if (ix[i] != i)
			continue;

		const MChain* chain = chains[i];
		const vector<MResidue*>& residues(chain->GetResidues());
		
		string& seq = seqset[i];
		assert(not seq.empty());
		seqlength += seq.length();

		vector<hit_ptr> c_hits;
		vector<res_ptr> c_res;

		ChainToHits(inDatabank, alignments[i], c_hits, residues, c_res);
		
		if (not c_hits.empty() and not c_res.empty())
		{
			assert(c_res.size() == residues.size());
			
			for (uint32 i = 0; i < c_res.size(); ++i)
			{
				assert(kResidueInfo[residues[i]->GetType()].code == c_res[i]->letter);
				
				if (i > 0 and residues[i - 1]->GetNumber() + 1 != residues[i]->GetNumber())
					++seqNr;
					
				c_res[i]->seqNr = seqNr;
				++seqNr;

				c_res[i]->chain = chain->GetChainID();
				c_res[i]->dssp = ResidueToDSSPLine(inProtein, *residues[i]).substr(5, 34);
			}
			
			msas.push_back(MSAInfo(seq, c_hits, c_res));
			hits.insert(hits.end(), c_hits.begin(), c_hits.end());
			res.insert(res.end(), c_res.begin(), c_res.end());
		}
	}

	sort(hits.begin(), hits.end(), compare_hit());
	if (hits.size() > 9999)
		hits.erase(hits.begin() + 9999, hits.end());
	
	uint32 nr = 1;
	foreach (hit_ptr h, hits)
		h->nr = nr++;
	
	string usedChains;
	for (int i = 0; i < ix.size(); ++i)
	{
		if (ix[i] == i)
		{
			if (not usedChains.empty())
				usedChains += ',';
			usedChains += chains[i]->GetChainID();
		}
	}
	
	CreateHSSPOutput(inProtein, inDatabank->GetVersion(), seqlength,
		chains.size(), msas.size(), usedChains, hits, res, msas, outHSSP);
}

}
