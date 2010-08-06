// i/o code for mas alignments and sequences

#include "mas.h"

#include <iostream>

#include <boost/filesystem/fstream.hpp>
#include <boost/foreach.hpp>
#define foreach BOOST_FOREACH
#include <boost/iostreams/device/array.hpp>
#include <boost/iostreams/stream.hpp>
#include <boost/format.hpp>

#include "ioseq.h"
#include "utils.h"

using namespace std;
namespace fs = boost::filesystem;
namespace io = boost::iostreams;

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

