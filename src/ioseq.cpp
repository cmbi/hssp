// i/o code for mas alignments and sequences

#include "mas.h"

#include <iostream>

#include <boost/filesystem/fstream.hpp>
#include <boost/foreach.hpp>
#define foreach BOOST_FOREACH
#include <boost/iostreams/device/array.hpp>
#include <boost/iostreams/stream.hpp>
#include <boost/format.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/regex.hpp>

#include "ioseq.h"
#include "utils.h"

using namespace std;
namespace fs = boost::filesystem;
namespace io = boost::iostreams;
namespace ba = boost::algorithm;

// --------------------------------------------------------------------

void readFasta(fs::path path, vector<entry>& seq)
{
	fs::ifstream file(path);
	if (not file.is_open())
		throw mas_exception(boost::format("input file '%1%' not opened") % path.string());
	
	string id, s;
	
	for (;;)
	{
		string line;
		getline(file, line);
		if (line.empty() or line[0] == '>')
		{
			if (line.empty() and not file.eof())
				continue;
			
			if (not (id.empty() or s.empty()))
			{
				if (VERBOSE)
					cout << "Sequence " << seq.size() + 1 << ": "
						 << id << string(20 - id.length(), ' ') << s.length() << " aa" << endl;
				seq.push_back(entry(seq.size(), id, encode(s)));
			}
			id.clear();
			s.clear();
			
			if (line.empty())
				break;
			
			string::size_type w = line.find(' ');
			if (w != string::npos)
				id = line.substr(1, w - 1);
			else
				id = line.substr(1);
			
			s.clear();
			continue;
		}
		
		s += line;
	}
}

void readAlignmentFromHsspFile(
	fs::path		path,
	char&			chainID,
	vector<entry>&	seq)
{
	seq.clear();

	fs::ifstream file(path);
	if (not file.is_open())
		throw mas_exception(boost::format("input file '%1%' not opened") % path.string());
	
	string line;
	
	// first line, should be something like "HSSP   bla bla bla"
	getline(file, line);
	if (not ba::starts_with(line, "HSSP"))
		throw mas_exception("file is not an HSSP file, does not start with HSSP");
	
	// second line, should contain PDBID
	getline(file, line);
	if (not ba::starts_with(line, "PDBID") or line.length() < 13)
		throw mas_exception("file is not a valid HSSP file, PDBID missing");
	
	entry pdb(0, line.substr(11), sequence());
	seq.push_back(pdb);
	
	string header;
	uint32 seqLength = 0, nchain = 0, nalign = 0, chain = 0;
	
	while (not file.eof())
	{
		getline(file, line);
		
		if (ba::starts_with(line, "## "))
			break;
		
		if (ba::starts_with(line, "HEADER     "))
			header = line.substr(11);
		else if (ba::starts_with(line, "SEQLENGTH  "))
			seqLength = boost::lexical_cast<uint32>(ba::trim_copy(line.substr(11, 4)));
		else if (ba::starts_with(line, "NCHAIN     "))
			nchain = boost::lexical_cast<uint32>(ba::trim_copy(line.substr(11, 4)));
		else if (ba::starts_with(line, "NALIGN     "))
			nalign = boost::lexical_cast<uint32>(ba::trim_copy(line.substr(11, 4)));
		else if (ba::starts_with(line, "KCHAIN     "))
		{
			vector<string> chains;
			string l = line.substr(48);
			ba::split(chains, l, ba::is_any_of(","));
			
			if (chains.empty())
				throw mas_exception("Invalid KCHAIN line in HSSP file");
			
			if (chainID == 0)
				chainID = chains.front()[0];
			else
			{
				for (chain = 0; chain < chains.size(); ++chain)
				{
					if (chainID == chains[chain][0])
						break;
				}
				
				if (chain == chains.size())
					throw mas_exception(boost::format("Chain %c not found in HSSP file") % chainID);
			}
		}
	}
	
	if (nalign < 1)
		throw mas_exception("invalid HSSP file, number of alignments missing");
	
	if (seqLength < 1)
		throw mas_exception("invalid HSSP file, sequence length missing");
	
	if (nchain < 1)
		throw mas_exception("invalid HSSP file, nchain missing");
	
	// second part, collect proteins
	
	if (line != "## PROTEINS : EMBL/SWISSPROT identifier and alignment statistics")
		throw mas_exception("invalid or unsupported HSSP file, ## PROTEINS line does match expected value");

	getline(file, line);
	if (not ba::starts_with(line, "  NR."))
		throw mas_exception("invalid HSSP file, expected line starting with NR.");
	
	seq.reserve(nalign);
	
	for (uint32 nr = 0; file.eof() == false and nr < nalign; ++nr)
	{
		getline(file, line);
		
		if (line.length() < 92)
			throw mas_exception("invalid HSSP file, protein line too short");
		
		entry e(nr + 1, line.substr(8, 12), sequence());
		ba::trim(e.m_id);
		
		string len = line.substr(74, 4);
		ba::trim(len);
		e.m_seq.reserve(boost::lexical_cast<uint32>(len));
		
		seq.push_back(e);
	}
	
	assert(seq.size() == nalign + 1);
	
	getline(file, line);
	
	uint32 alignment = 1;
	while (file.eof() == false and alignment < nalign)
	{
		if (not ba::starts_with(line, "## ALIGNMENTS"))
			throw mas_exception("invalid HSSP file, missing ## ALIGNMENTS line");

		string a_from = line.substr(14, 4);	ba::trim(a_from);
		string a_to = line.substr(21, 4); ba::trim(a_to);
		
		uint32 a_f = boost::lexical_cast<uint32>(a_from);
		uint32 a_t = boost::lexical_cast<uint32>(a_to);
		
		if (a_f != alignment or a_t < a_f or a_t > nalign)
			throw mas_exception("Invalid HSSP file, incorrect number of alignments");
		
		alignment = a_t + 1;
		
		getline(file, line);	// SeqNo line
		
		string pdb_seq;
		vector<string> s(a_t - a_f + 1);
		vector<uint16> pdbnrs;
		
		do {
			getline(file, line);
			
			int32 n = static_cast<int32>(line.length()) - 51;
			
			if (chainID == 0)			// store chain ID
				chainID = line[12];

			if (line[12] == chainID)
			{
				string pdbno = line.substr(7, 4);
				ba::trim(pdbno);
				uint16 pdbno_value = boost::lexical_cast<uint16>(pdbno);
				
				pdbnrs.push_back(pdbno_value);
				
				pdb_seq += line[14];
				
				for (int32 i = 0; i < n; ++i)
				{
					char a = line[51 + i];
					if (a == ' ' or a == '.')
						a = '-';

					s[i] += a;
					seq[a_f + i].m_positions.push_back(pdbno_value);
				}
			}
		}
		while (not file.eof() and not (ba::starts_with(line, "##") or ba::starts_with(line, "//")));
		
		for (uint32 i = 0; i < a_t - a_f + 1; ++i)
		{
			seq[a_f + i].m_seq = encode(s[i]);
			assert(seq[a_f + i].m_seq.length() == seq[a_f + i].m_positions.size());
		}
		
		if (seq.front().m_seq.empty())
		{
			seq.front().m_seq = encode(pdb_seq);
			seq.front().m_positions = pdbnrs;
		}
		else if (decode(seq.front().m_seq) != pdb_seq)
		{
			cout << endl << decode(seq.front().m_seq) << endl
				 << pdb_seq << endl
				 << endl; 
			
			throw mas_exception("Invalid HSSP file, inconsistent PDB sequence");
		}
	}
	
	// process ## INSERTION LIST, if any
	
	while (not file.eof() and ba::starts_with(line, "## ") and
		not ba::starts_with(line, "## INSERTION LIST"))
	{
		do
		{
			getline(file, line);
		}
		while (not (file.eof() or ba::starts_with(line, "## ")));
	}
	
	if (ba::starts_with(line, "## INSERTION LIST"))
	{
		getline(file, line);
		
		for (;;)
		{
			getline(file, line);
			
			if (file.eof() or ba::starts_with(line, "//"))
				break;

			uint32 l = boost::lexical_cast<uint32>(ba::trim_copy(line.substr(20, 4)));
			
			if (line.length() != 25 + l + 2)
				throw mas_exception("Invalid HSSP file, incorrect insertion length");
			
			sequence ins = encode(line.substr(26, l));
			
			uint32 p = boost::lexical_cast<uint32>(ba::trim_copy(line.substr(8, 4))) - 1;
			uint32 a = boost::lexical_cast<uint32>(ba::trim_copy(line.substr(2, 4)));
			
			seq[a].m_seq.insert(p, ins);
			seq[a].m_positions.insert(seq[a].m_positions.begin() + p, l, 0);
		}
	}
	
	seq.erase(
		remove_if(seq.begin(), seq.end(), boost::bind(&entry::length, _1) == 0),
		seq.end());
	
//	foreach (const entry& e, seq)
//	{
//		cout << e.m_id << endl
//			 << decode(e.m_seq) << endl;
//		
//		foreach (uint16 p, e.m_positions)
//			cout << p % 10;
//
//		cout << endl << endl;
//	}
	
}

void readWhatifMappingFile(fs::path path, vector<entry>& seq)
{
	seq.clear();

	fs::ifstream file(path);
	if (not file.is_open())
		throw mas_exception(boost::format("input file '%1%' not opened") % path.string());
	
	string line;
	
	while (not file.eof())
	{
		getline(file, line);
		if (not ba::starts_with(line, "Sequence name: "))
			continue;
		break;
	}
	
	while (ba::starts_with(line, "Sequence name: "))
	{
		string id = line.substr(15);
		string s;
		vector<uint16> pos;
		
		// skip the description line
		getline(file, line);
		
		while (not (file.eof() or ba::starts_with(line, "Sequence name: ")))
		{
			getline(file, line);

			if (line.length() != 16 or line[4] != ' ' or line[6] != ' ' or line[11] != ' ')
				continue;
			
			s += line[5];
			
			string nr = line.substr(12);
			ba::trim(nr);
			
			if (nr == "----")
				pos.push_back(0);
			else
				pos.push_back(boost::lexical_cast<uint16>(nr));
		}
		
		if (not s.empty())
		{
			entry e(seq.size(), id, encode(s));
			e.m_positions = pos;
			seq.push_back(e);
		}
	}
}

void readFamilyIdsFile(fs::path path, vector<entry>& seq)
{
	seq.clear();

	fs::ifstream file(path);
	if (not file.is_open())
		throw mas_exception(boost::format("input file '%1%' not opened") % path.string());
	
	fs::path dir = path.parent_path();
	
	while (not file.eof())
	{
		string id;
		getline(file, id);
		
		if (id.empty())
			continue;
		
		fs::ifstream data(dir / (id + ".mapping"));
		if (not data.is_open())
			throw mas_exception(boost::format("Failed to open mapping file for protein %1%") % id);
		
		string line, s;
		vector<uint16> pos;

		while (not data.eof())
		{
			getline(data, line);
	
			if (line.length() < 3 or line[1] != '\t')
				continue;
			
			s += line[0];
			pos.push_back(boost::lexical_cast<uint16>(line.substr(2)));
		}
		
		if (not s.empty())
		{
			entry e(seq.size(), id, encode(s));
			e.m_positions = pos;
			seq.push_back(e);
		}
	}
}

// --------------------------------------------------------------------

void report_in_fasta(const vector<entry*>& alignment, ostream& os)
{
	foreach (const entry* e, alignment)
	{
		os << '>' << e->m_id << endl;
		
		uint32 o = 0;
		while (o < e->m_seq.length())
		{
			uint32 n = e->m_seq.length() - o;
			if (n > 72)
				n = 72;
			
			os << decode(e->m_seq.substr(o, n)) << endl;
			o += n;
		}
	}
}

void report_in_clustalw(const vector<entry*>& alignment, ostream& os)
{
	os << "CLUSTAL FORMAT for MaartensAlignment" << endl;

	uint32 nseq = alignment.size();
	uint32 len = alignment[0]->m_seq.length();
	uint32 offset = 0;
	while (offset < len)
	{
		uint32 n = alignment[0]->m_seq.length() - offset;
		if (n > 60)
			n = 60;
		
		struct {
			uint32		cnt[20];
		} dist[60] = {};
		
		foreach (const entry* e, alignment)
		{
			sequence ss = e->m_seq.substr(offset, n);
			
			for (uint32 i = 0; i < n; ++i)
			{
				aa ri = ss[i];
				if (ri < 20)
					dist[i].cnt[ri] += 1;
			}

			string id = e->m_id;
			if (id.length() > 15)
				id = id.substr(0, 12) + "...";
			else if (id.length() < 15)
				id += string(15 - id.length(), ' ');
			
			os << id << ' ' << decode(ss) << endl;
		}
		
		string scores(n, ' ');
		for (uint32 i = 0; i < n; ++i)
		{
			uint32 strong[9] = {};
			const char* kStrongGroups[9] = {
				"STA", "NEQK", "NHQK", "NDEQ", "QHRK", "MILV", "MILF", "HY", "FYW"
			};
			
			uint32 weak[11] = {};
			const char* kWeakGroups[11] = {
				"CSA", "ATV", "SAG", "STNK", "STPA", "SGND", "SNDEQK",
				"NDEQHK", "NEQHRK", "FVLIM", "HFY"
			};
			
			for (uint32 r = 0; r < 20; ++r)
			{
				if (dist[i].cnt[r] == alignment.size())
				{
					scores[i] = '*';
					break;
				}
				
				for (uint32 g = 0; g < 9; ++g)
				{
					if (strchr(kStrongGroups[g], kAA[r]) != NULL)
						strong[g] += dist[i].cnt[r];
				}

				for (uint32 g = 0; g < 11; ++g)
				{
					if (strchr(kWeakGroups[g], kAA[r]) != NULL)
						weak[g] += dist[i].cnt[r];
				}
			}
			
			for (uint32 g = 0; scores[i] == ' ' and g < 9; ++g)
			{
				if (strong[g] == alignment.size())
					scores[i] = ':';
			}

			for (uint32 g = 0; scores[i] == ' ' and g < 11; ++g)
			{
				if (weak[g] == alignment.size())
					scores[i] = '.';
			}
		}
		
		os << string(16, ' ') << scores << endl;

		if (not alignment.front()->m_positions.empty())
		{
			string pos_nrs(n, ' ');
			for (uint32 i = 0; i < n; ++i)
			{
				if (alignment.front()->m_positions[offset + i] != 0)
					pos_nrs[i] = '!';
			}
			
			os << string(16, ' ') << pos_nrs << endl;
		}
		
		offset += n;
		os << endl;
	}
}

void report(const vector<entry*>& alignment, ostream& os, const string& format)
{
	if (format == "fasta")
		report_in_fasta(alignment, os);
	else if (format == "clustalw")
		report_in_clustalw(alignment, os);
	else
		throw mas_exception(boost::format("Unknown output format %1%") % format);
}

