// align.cpp - simple attempt to write a multiple sequence alignment application
//

#include "align.h"

#include <iostream>
#include <iomanip>
#include <string>
#include <sstream>
#include <limits>
#include <cmath>
#include <numeric>

#include <boost/program_options.hpp>
#include <boost/format.hpp>
#include <boost/filesystem/fstream.hpp>
#include <boost/foreach.hpp>
#include <boost/tr1/tuple.hpp>
#define foreach BOOST_FOREACH
#include <boost/iostreams/device/array.hpp>
#include <boost/iostreams/stream.hpp>
#include <boost/bind.hpp>
#include <boost/algorithm/string.hpp>

#include "matrix.h"

using namespace std;
using namespace tr1;
namespace fs = boost::filesystem;
namespace po = boost::program_options;
namespace io = boost::iostreams;
namespace ba = boost::algorithm;

int DEBUG = 0, VERBOSE = 0;

// --------------------------------------------------------------------

aa kAA_Reverse[256];

struct entry
{
					entry(uint32 nr, const string& id, const sequence& seq, float weight = 1.0f)
						: m_nr(nr)
						, m_id(id)
						, m_seq(seq)
						, m_weight(weight) {}

	uint32			nr() const						{ return m_nr; }
	float			weight() const					{ return m_weight; }

	uint32			m_nr;
	string			m_id;
	sequence		m_seq;
	float			m_weight;
};

// ah, too bad those lamda's are not supported yet...
struct sum_weight
{
	float operator()(float sum, const entry* e) const { return sum + e->m_weight; }
};

// --------------------------------------------------------------------

my_bad::my_bad(const string& msg)
{
	snprintf(m_msg, sizeof(m_msg), "%s", msg.c_str());
}

my_bad::my_bad(const boost::format& msg)
{
	snprintf(m_msg, sizeof(m_msg), "%s", msg.str().c_str());
}

// --------------------------------------------------------------------

template<typename T>
class distance_matrix
{
  public:
	typedef T value_type;
	
					distance_matrix(uint32 n);
	virtual			~distance_matrix();
	
	value_type		operator()(uint32 i, uint32 j) const;
	value_type&		operator()(uint32 i, uint32 j);
	
	// erase two rows, add one at the end (for neighbour joining)
	void			erase_2(uint32 i, uint32 j);

	void			print(ostream& os) const;

  private:
	value_type*		m_data;
	uint32			m_n;
};

template<typename T>
distance_matrix<T>::distance_matrix(uint32 n)
	: m_n(n)
{
	m_data = new value_type[(m_n * (m_n - 1)) / 2];
}

template<typename T>
distance_matrix<T>::~distance_matrix()
{
	delete[] m_data;
}

template<typename T>
inline
T& distance_matrix<T>::operator()(uint32 i, uint32 j)
{
	if (i > j)
		swap(i, j);
	
	assert(j < m_n); assert(i != j);
	return m_data[(j * (j - 1)) / 2 + i];
}

template<typename T>
inline
T distance_matrix<T>::operator()(uint32 i, uint32 j) const
{
	if (i > j)
		swap(i, j);
	
	assert(j < m_n); assert(i != j);
	return m_data[(j * (j - 1)) / 2 + i];
}

template<typename T>
void distance_matrix<T>::erase_2(uint32 di, uint32 dj)
{
	uint32 s = 0, d = 0;
	for (uint32 i = 0; i < m_n; ++i)
	{
		for (uint32 j = 0; j < i; ++j)
		{
			if (i != di and j != dj and i != dj and j != di)
			{
				if (s != d)
					m_data[d] = m_data[s];
				++d;
			}

			++s;
		}
	}
	
	--m_n;
}

template<typename T>
void distance_matrix<T>::print(ostream& os) const
{
	for (uint32 y = 1; y < m_n; ++y)
	{
		cout << setw(5) << y;
		
		for (uint32 x = 0; x < y; ++x)
			cout << (boost::format("  %1.2f") % operator()(x, y));
		
		cout << endl;
	}
}

template<typename T>
ostream& operator<<(ostream& lhs, const distance_matrix<T>& rhs)
{
	rhs.print(lhs);
	return lhs;
}

// --------------------------------------------------------------------

struct base_node
{
	virtual				~base_node() {}

	virtual void		print(ostream& s) = 0;	

	virtual base_node*	left() const		{ return 0; }
	virtual base_node*	right() const		{ return 0; }

	virtual void		add_weight(float w) = 0;
	virtual uint32		leaf_count() const	{ return 1; }
};

ostream& operator<<(ostream& lhs, base_node& rhs)
{
	rhs.print(lhs);
	return lhs;
}

struct joined_node : public base_node
{
						joined_node(base_node* left, base_node* right,
							float d_left, float d_right);

	virtual				~joined_node()
						{
							delete m_left;
							delete m_right;
						}

	virtual void		print(ostream& s)
						{
							s << '(' << endl
							  << *m_left << (boost::format(":%1.4f") % m_d_left) << ',' << endl
							  << *m_right << (boost::format(":%1.4f") % m_d_right) << ')' << endl;
						}

	virtual base_node*	left() const		{ return m_left; }
	virtual base_node*	right() const		{ return m_right; }

	virtual void		add_weight(float w);
	virtual uint32		leaf_count() const	{ return m_leaf_count; }
	
	base_node*			m_left;
	base_node*			m_right;
	float				m_d_left;
	float				m_d_right;
	uint32				m_leaf_count;
};

struct leaf_node : public base_node
{
						leaf_node(entry& e)
							: m_entry(e)
						{
							m_entry.m_weight = 0;
						}

	virtual void		print(ostream& s)
						{
							s << m_entry.m_id;
							if (VERBOSE)
								s << '(' << m_entry.m_weight << ')';
						}

	virtual void		add_weight(float w)
						{
							m_entry.m_weight += w;
						}

	entry&				m_entry;
};

joined_node::joined_node(base_node* left, base_node* right, float d_left, float d_right)
	: m_left(left)
	, m_right(right)
	, m_d_left(d_left)
	, m_d_right(d_right)
	, m_leaf_count(left->leaf_count() + right->leaf_count())
{
	m_left->add_weight(d_left / m_left->leaf_count());
	m_right->add_weight(d_right / m_right->leaf_count());
}

void joined_node::add_weight(float w)
{
	m_left->add_weight(w);
	m_right->add_weight(w);
}

// --------------------------------------------------------------------

void report(const vector<entry*>& alignment)
{
	cout << "CLUSTAL FORMAT for MaartensAlignment" << endl;

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
			
			cout << id << ' ' << decode(ss) << endl;
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
		
		cout << string(16, ' ') << scores << endl;
		
		offset += n;
		cout << endl;
	}
}

// --------------------------------------------------------------------
// distance is calculated as 1 minus the fraction of identical residues

float calculateDistance(const entry& a, const entry& b)
{
	int32 x, dimX = a.m_seq.length();
	int32 y, dimY = b.m_seq.length();
	
	matrix<float>	B(dimX, dimY);
	matrix<float>	Ix(dimX, dimY);
	matrix<float>	Iy(dimX, dimY);
	matrix<uint16>	id(dimX, dimY);
	
	Ix(0, 0) = 0;
	Iy(0, 0) = 0;
	
	static const substitution_matrix smat("GONNET250");
	static const float gapOpen = 10;
	static const float gapExtend = 0.2;
	
	float high = numeric_limits<float>::min();
	uint32 highX, highY, highId;

	for (x = 0; x < dimX; ++x)
	{
		for (y = 0; y < dimY; ++y)
		{
			float Ix1 = 0; if (x > 0) Ix1 = Ix(x - 1, y);
			float Iy1 = 0; if (y > 0) Iy1 = Iy(x, y - 1);

			// (1)
			float M = smat(a.m_seq[x], b.m_seq[y]);
			if (x > 0 and y > 0)
				M += B(x - 1, y - 1);

			float s;
			uint16 i = 0;
			if (a.m_seq[x] == b.m_seq[y])
				i = 1;

			if (M >= Ix1 and M >= Iy1)
			{
				if (x > 0 and y > 0)
					i += id(x - 1, y - 1);
				s = M;
			}
			else if (Ix1 >= Iy1)
			{
				if (x > 0)
					i += id(x - 1, y);
				s = Ix1;
			}
			else
			{
				if (y > 0)
					i += id(x, y - 1);
				s = Iy1;
			}
			
			B(x, y) = s;
			id(x, y) = i;
			
			if ((x == dimX - 1 or y == dimY - 1) and high < s)
			{
				high = s;
				highX = x;
				highY = y;
				highId = i;
			}

			// (3)
			Ix(x, y) = max(M - gapOpen, Ix1 - gapExtend);
			
			// (4)
			Iy(x, y) = max(M - gapOpen, Iy1 - gapExtend);
		}
	}
	
	float result = 1.0f - float(highId) / max(dimX, dimY);
	
	if (VERBOSE)
		cout << (boost::format("Sequences (%1$d:%2$d) Aligned. Score: %3$4.2f") % (a.m_nr + 1) % (b.m_nr + 1) % result) << endl;
	
	return result;
}

// --------------------------------------------------------------------

string decode(const sequence& s)
{
	string result;
	result.reserve(s.length());
	
	foreach (aa a, s)
		result.push_back(kAA[a]);

	return result;
}

sequence encode(const string& s)
{
	static bool sInited = false;
	static uint8 kAA_Reverse[256];
	
	if (not sInited)
	{
		// init global reverse mapping
		for (uint32 a = 0; a < 256; ++a)
			kAA_Reverse[a] = 255;
		for (uint32 a = 0; a < sizeof(kAA); ++a)
			kAA_Reverse[kAA[a]] = a;
	}
	
	sequence result;
	result.reserve(s.length());

	foreach (char r, s)
	{
		aa rc = kAA_Reverse[static_cast<uint8>(r)];
		if (rc < sizeof(kAA))
			result.push_back(rc);
	}
	
	return result;
}

void readFasta(fs::path path, vector<entry>& seq)
{
	fs::ifstream file(path);
	if (not file.is_open())
		throw my_bad(boost::format("input file '%1%' not opened") % path.string());
	
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

void joinNeighbours(distance_matrix<float>& d, vector<base_node*>& tree)
{
	uint32 r = tree.size();
	
	while (r > 2)
	{
		// calculate the sums first
		vector<float> sum(r);
		for (uint32 i = 1; i < r; ++i)
		{
			for (uint32 j = 0; j < i; ++j)
			{
				float dij = d(i, j);
				sum[i] += dij;
				sum[j] += dij;
			}
		}
		
		// calculate Q, or in fact, the position of the minimum in Q
		uint32 min_i = 0, min_j = 0;
		float m = numeric_limits<float>::max();
		
		for (uint32 i = 1; i < r; ++i)
		{
			for (uint32 j = 0; j < i; ++j)
			{
				float v = d(i, j) - (sum[i] + sum[j]) / (r - 2);
	
				if (m > v)
				{
					min_i = i;
					min_j = j;
					m = v;
				}
			}
		}
		
		// distance to joined node
		float d_i, d_j;
		float half_dij = d(min_i, min_j) / 2;
		
		d_i = half_dij + abs(sum[min_i] - sum[min_j]) / (2 * (r - 2));
		d_j = d(min_i, min_j) - d_i;
	
		if (d_i > d_j and tree[min_i]->leaf_count() > tree[min_j]->leaf_count())
			swap(d_i, d_j);
	
		joined_node* jn = new joined_node(tree[min_i], tree[min_j], d_i, d_j);
		assert(min_j < min_i);
		tree.erase(tree.begin() + min_i);
		tree.erase(tree.begin() + min_j);
		tree.push_back(jn);
		
		// distances to other nodes
		vector<float> dn; dn.reserve(r - 2);
		for (uint32 x = 0; x < r; ++x)
		{
			if (x == min_i or x == min_j)
				continue;
			dn.push_back((abs(d(x, min_i) - d_i) + abs(d(x, min_j) - d_j)) / 2);
		}
		
		// fill new distance matrix
		d.erase_2(min_i, min_j);
		--r;
		for (uint32 x = 0; x < r - 1; ++x)
			d(x, r - 1) = dn[x];
	}

	assert(r == 2); assert(tree.size() == 2);

	joined_node* root = new joined_node(tree[0], tree[1], d(0, 1) / 2, d(0, 1) / 2);
	tree.clear();
	tree.push_back(root);
}

// --------------------------------------------------------------------

class GuideTreeParser
{
  public:
						GuideTreeParser(istream& data, map<string,leaf_node*>& m)
							: m_data(data)
							, m_map(m)
						{
							getNextToken();
						}
	
	base_node*			parse();

  private:

	enum GuideTreeToken {
		gtt_Undef,
		gtt_Open,
		gtt_Close,
		gtt_End,
		gtt_ID,
		gtt_Colon,
		gtt_Comma,
		gtt_Weight
	};
	
	void				getNextToken();
	void				match(GuideTreeToken token);
	
	tr1::tuple<base_node*,float>
						parseNode();
	base_node*			parseGroup();
	
	istream&			m_data;
	map<string,leaf_node*>&
						m_map;
	float				m_value;
	string				m_token;
	GuideTreeToken		m_lookahead;
};

void GuideTreeParser::getNextToken()
{
	m_lookahead = gtt_Undef;
	m_token.clear();
	m_value = 0;
	
	enum State {
		st_Start,
		st_ID,
		st_Number,
		st_Fraction
	} state = st_Start, start = st_Start;
	
	while (m_lookahead == gtt_Undef)
	{
		char ch = 0;
		m_data.get(ch);
		
		m_token += ch;
		
		switch (state)
		{
			case st_Start:
				switch (ch)
				{
					case ' ':
					case '\r':
					case '\n':
					case '\t':	m_token.clear();			break;
					case '(':	m_lookahead = gtt_Open;		break;
					case ')':	m_lookahead = gtt_Close;	break;
					case ':':	m_lookahead = gtt_Colon;	break;
					case ',':	m_lookahead = gtt_Comma;	break;
					case ';':	m_lookahead = gtt_End;		break;
					default:
						if (isdigit(ch) or ch == '-')
							state = st_Number;
						else if (isalnum(ch) or ch == '_')
							state = st_ID;
						else
							throw my_bad(boost::format("unexpected character '%1%' in guide tree") % ch);
						break;
				}
				break;
			
			case st_Number:
				if (ch == '.')
					state = st_Fraction;
				else if (not isdigit(ch))
				{
					m_data.unget();
					m_token.erase(m_token.end() - 1);
					m_lookahead = gtt_Weight;
				}
				break;
			
			case st_Fraction:
				if (not isdigit(ch))
				{
					m_data.unget();
					m_token.erase(m_token.end() - 1);
					m_lookahead = gtt_Weight;
				}
				break;
			
			case st_ID:
				if (not isalnum(ch) and ch != '_')
				{
					m_data.unget();
					m_token.erase(m_token.end() - 1);
					m_lookahead = gtt_ID;
				}
				break;
		}
	}
}

void GuideTreeParser::match(GuideTreeToken token)
{
	if (token == m_lookahead)
		getNextToken();
	else
		throw my_bad(boost::format("invalid guide tree, expected %1% but found %2% ('%3%')") %
			int(token) % int(m_lookahead) % m_token);
}

tr1::tuple<base_node*,float> GuideTreeParser::parseNode()
{
	base_node* n = NULL;
	float w = 0;
	
	if (m_lookahead == gtt_Open)
		n = parseGroup();
	else
	{
		n = m_map[m_token];
		if (n == NULL)
			throw my_bad(boost::format("guide tree contains unknown id %1%") % m_token);
		match(gtt_ID);
	}
	
	if (m_lookahead == gtt_Colon)
	{
		match(gtt_Colon);
		w = boost::lexical_cast<float>(m_token);
		match(gtt_Weight);
	}
	
	return tr1::make_tuple(n, w);
}

base_node* GuideTreeParser::parseGroup()
{
	base_node* na = NULL;
	base_node* nb = NULL;
	float wa = 0, wb = 0;
	
	match(gtt_Open);
	
	tr1::tie(na, wa) = parseNode();
	
	while (m_lookahead == gtt_Comma)
	{
		match(gtt_Comma);
		
		tr1::tie(nb, wb) = parseNode();
		
		na = new joined_node(na, nb, wa, wb);
		wa = wb;
	}
	
	match(gtt_Close);
	
	return na;
}

base_node* GuideTreeParser::parse()
{
	base_node* result = NULL;
	if (m_lookahead == gtt_Open)
		result = parseGroup();
	if (m_lookahead != gtt_End)
		throw my_bad("invalid guide tree file, missing semicolon at end");
	return result;
}

void useGuideTree(const string& guide, vector<base_node*>& tree)
{
	uint32 r = tree.size();
	
	map<string,leaf_node*> m;
	foreach (base_node* n, tree)
	{
		leaf_node* leaf = static_cast<leaf_node*>(n);
		m[leaf->m_entry.m_id] = leaf;
	}
	
	fs::ifstream file(guide);
	if (not file.is_open())
		throw my_bad("failed to open guide tree");
	
	tree.clear();
	GuideTreeParser parser(file, m);
	tree.push_back(parser.parse());
}

// --------------------------------------------------------------------

float score(const vector<entry*>& a, const vector<entry*>& b,
	uint32 ix_a, uint32 ix_b, const substitution_matrix& mat)
{
	float result = 0;
	
	foreach (const entry* ea, a)
	{
		foreach (const entry* eb, b)
		{
			assert(ix_a < ea->m_seq.length());
			assert(ix_b < eb->m_seq.length());
			
			aa ra = ea->m_seq[ix_a];
			aa rb = eb->m_seq[ix_b];
			
			if (ra != kSignalGapCode and rb != kSignalGapCode)
				result += ea->m_weight * eb->m_weight * mat(ra, rb);
		}
	}
	
	return result / (a.size() * b.size());
}

float percentIdentity(const sequence& a, const sequence& b)
{
	float result = 0;
	uint32 count = 0, total;
	
	total = min(a.length(), b.length());
	
	for (uint32 ix = 0; ix < total; ++ix)
	{
		if (a[ix] == b[ix])
			++count;
	}
	
	if (total > 0)
		result = 100.0f * count / total;
	
	return result;
}

// don't ask me, but looking at the clustal code, they substract 0.2 from the table
// as mentioned in the article in NAR.
const float kResidueSpecificPenalty[20] = {
	1.13 - 0.2,		// A
	0.72 - 0.2,		// R
	0.63 - 0.2,		// N
	0.96 - 0.2,		// D
	1.13 - 0.2,		// C
	1.07 - 0.2,		// Q
	1.31 - 0.2,		// E
	0.61 - 0.2,		// G
	1.00 - 0.2,		// H
	1.32 - 0.2,		// I
	1.21 - 0.2,		// L
	0.96 - 0.2,		// K
	1.29 - 0.2,		// M
	1.20 - 0.2,		// F
	0.74 - 0.2,		// P
	0.76 - 0.2,		// S
	0.89 - 0.2,		// T
	1.23 - 0.2,		// W
	1.00 - 0.2,		// Y
	1.25 - 0.2		// V
};

void adjust_gp(vector<float>& gop, vector<float>& gep, const vector<entry*>& seq)
{
	assert(gop.size() == seq.front()->m_seq.length());

	vector<uint32> gaps(gop.size());
	vector<bool> hydrophilic_stretch(gop.size(), false);
	vector<float> residue_specific_penalty(gop.size());

	foreach (const entry* e, seq)
	{
		const sequence& s = e->m_seq;
		
		for (uint32 ix = 0; ix < gop.size(); ++ix)
		{
			aa r = s[ix];
			
			if (r == kSignalGapCode)
				gaps[ix] += 1;

			// residue specific gap penalty
			if (r < kAACount)
				residue_specific_penalty[ix] += kResidueSpecificPenalty[r];
			else
				residue_specific_penalty[ix] += 1.0f;
		}
		
		// find a run of 5 hydrophilic residues
		static const boost::function<bool(aa)> is_hydrophilic = ba::is_any_of(encode("DEGKNQPRS"));
		
		for (uint32 si = 0, i = 0; i <= gop.size(); ++i)
		{
			if (i == gop.size() or is_hydrophilic(s[i]) == false)
			{
				if (i >= si + 5)
				{
					for (uint32 j = si; j < i; ++j)
						hydrophilic_stretch[j] = true;
				}
				si = i + 1;
			}
		}
	}
	
	for (int32 ix = 0; ix < static_cast<int32>(gop.size()); ++ix)
	{
		// if there is a gap, lower gap open cost
		if (gaps[ix] > 0)
		{
			gop[ix] *= 0.3f * ((seq.size() - gaps[ix]) / float(seq.size()));
			gep[ix] /= 2;
		}
		
		// else if there is a gap within 8 residues, increase gap cost
		else
		{
			for (int32 d = 0; d < 8; ++d)
			{
				if (ix + d >= gaps.size() or gaps[ix + d] > 0 or
					ix - d < 0 or gaps[ix - d] > 0)
				{
					gop[ix] *= (2 + ((8 - d) * 2)) / 8.f;
					break;
				}
			}
			
			if (hydrophilic_stretch[ix])
				gop[ix] /= 3;
			else
				gop[ix] *= (residue_specific_penalty[ix] / seq.size());
		}
	}

	if (DEBUG > 2)
	{
		foreach (const entry* e, seq)
			cout << e->m_id << " (" << gop.size() << "; " << e->m_weight << ")" << endl;
		copy(gop.begin(), gop.end(), ostream_iterator<float>(cout, ";")); cout << endl;
	}
}

void align(
	const joined_node* node,
	vector<entry*>& a, vector<entry*>& b, vector<entry*>& c,
	const substitution_matrix_family& mat_fam,
	float gop, float gep)
{
	int32 x, dimX = a.front()->m_seq.length();
	int32 y, dimY = b.front()->m_seq.length();
	
	matrix<float> B(dimX, dimY);
	matrix<float> Ix(dimX, dimY);
	matrix<float> Iy(dimX, dimY);
	matrix<int8> tb(dimX, dimY);
	
	Ix(0, 0) = 0;
	Iy(0, 0) = 0;

	const substitution_matrix& smat = mat_fam(abs(node->m_d_left + node->m_d_right), true);

	uint32 minLength = dimX, maxLength = dimY;
	if (minLength > maxLength)
		swap(minLength, maxLength);
	
	float logmin = 1.0 / log10(minLength);
	float logdiff = 1.0 + 0.5 * log10(float(minLength) / maxLength);
	
	// initial gap open cost, 0.05f is the remaining magical number here...
	gop = (gop / (logdiff * logmin)) * smat.mismatch_average() * smat.scale_factor() * 0.05f;

	float avg_weight_a = accumulate(a.begin(), a.end(), 0.f, sum_weight()) / a.size();
	float avg_weight_b = accumulate(b.begin(), b.end(), 0.f, sum_weight()) / b.size();

	// position specific gap open costs
	// initial gap extend cost is adjusted for difference in sequence lengths
	vector<float> gop_a(dimX, gop * avg_weight_a),
		gep_a(dimX, gep * (1 + log10(float(dimX) / dimY)) * avg_weight_a);
	adjust_gp(gop_a, gep_a, a);
	
	vector<float> gop_b(dimY, gop * avg_weight_b),
		gep_b(dimY, gep * (1 + log10(float(dimY) / dimX)) * avg_weight_b);
	adjust_gp(gop_b, gep_b, b);
	
	float high = numeric_limits<float>::min();
	int32 highX, highY;
	
	for (x = 0; x < dimX; ++x)
	{
		for (y = 0; y < dimY; ++y)
		{
			float Ix1 = 0; if (x > 0) Ix1 = Ix(x - 1, y);
			float Iy1 = 0; if (y > 0) Iy1 = Iy(x, y - 1);
			
			// (1)
			float M = score(a, b, x, y, smat);
			if (x > 0 and y > 0)
				M += B(x - 1, y - 1);

			float s;
			if (M >= Ix1 and M >= Iy1)
			{
				tb(x, y) = 0;
				B(x, y) = s = M;
			}
			else if (Ix1 >= Iy1)
			{
				tb(x, y) = 1;
				B(x, y) = s = Ix1;
			}
			else
			{
				tb(x, y) = -1;
				B(x, y) = s = Iy1;
			}
			
			if ((x == dimX - 1 or y == dimY - 1) and high <= s)
			{
				high = s;
				highX = x;
				highY = y;
			}
			
			if (DEBUG > 8)
			{
				cout << "x: " << x << "; y: " << y << "; m: " << M << "; Ix1: " << Ix1 << "; Iy1: " << Iy1
					 << "; B=> " << B(x, y) << "; tb=> " << int(tb(x, y)) << endl;
			}

			// (3)
			Ix(x, y) = max(M - gop_a[x], Ix1 - gep_a[x]);
			
			// (4)
			Iy(x, y) = max(M - gop_b[y], Iy1 - gep_b[y]);
		}
	}
	
	// build the final alignment
	x = dimX - 1;
	y = dimY - 1;

	// first append end-gaps
	while (x > highX)
	{
		foreach (entry* e, b)
			e->m_seq += kSignalGapCode;
		--x;
	}

	while (y > highY)
	{
		foreach (entry* e, a)
			e->m_seq += kSignalGapCode;
		--y;
	}

	while (x >= 0 or y >= 0)
	{
		if (x == -1)
		{
			foreach (entry* e, a)
				e->m_seq.insert(e->m_seq.begin(), kSignalGapCode);
			--y;
		}
		else if (y >= 0 and tb(x, y) < 0)
		{
			foreach (entry* e, a)
				e->m_seq.insert(e->m_seq.begin() + x, kSignalGapCode);
			--y;
		}
		else if (y == -1)
		{
			foreach (entry* e, b)
				e->m_seq.insert(e->m_seq.begin(), kSignalGapCode);
			--x;
		}
		else if (x >= 0 and tb(x, y) > 0)
		{
			foreach (entry* e, b)
				e->m_seq.insert(e->m_seq.begin() + y, kSignalGapCode);
			--x;
		}
		else
		{
			assert(x >= 0); assert(y >= 0);
			--x;
			--y;
		}
	}
	
	copy(a.begin(), a.end(), back_inserter(c));
	copy(b.begin(), b.end(), back_inserter(c));
	
	if (DEBUG == 2)
		report(c);
	
	assert(a.front()->m_seq.length() == b.front()->m_seq.length());

	if (DEBUG > 7)
	{
		foreach (entry* e, c)
			cout << e->m_id << ": " << decode(e->m_seq) << endl;
		cout << endl;
		
		cout << "      ";
		for (int32 x = 0; x < dimX; ++x)
			cout << setw(8) << x << "   ";
		cout << endl;
		
		for (int32 y = 0; y < dimY; ++y)
		{
			cout << setw(6) << y;
			for (int32 x = 0; x < dimX; ++x)
				cout << ' ' << setw(7) << B(x, y) << ' ' << setw(2) << int(tb(x, y));
			cout << endl;
		}
		
		cout << endl;
	}
}

void createAlignment(base_node* node, vector<entry*>& alignment,
	const substitution_matrix_family& mat, float gop, float gep)
{
	if (node->left() and node->right())
	{
		vector<entry*> a, b;	
		createAlignment(node->left(), a, mat, gop, gep);
		createAlignment(node->right(), b, mat, gop, gep);
		
		align(static_cast<joined_node*>(node), a, b,
			alignment, mat, gop, gep);
		
		if (DEBUG)
			report(alignment);
	}
	else
	{
		leaf_node* n = dynamic_cast<leaf_node*>(node);
		if (n == 0)
			throw my_bad("internal error (not a leaf node)");
		
		alignment.push_back(&n->m_entry);
	}
}

void test()
{
	vector<entry*> s;
	
	s.push_back(new entry(1, "1", encode("HLTPEEKSAVTALWGKVN--VDEVGGEALGRLLVVYPWTQRFFESFGDL")));
	s.push_back(new entry(2, "2", encode("QLSGEEKAAVLALWDKVN--EEEVGGEALGRLLVVYPWTQRFFDSFGDL")));
	s.push_back(new entry(3, "3", encode("VLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLS")));
	s.push_back(new entry(4, "4", encode("VLSAADKTNVKAAWSKVGGHAGEYGAEALERMFLGFPTTKTYFPHFDLS")));
	
	uint32 n = s.front()->m_seq.length();
	
	vector<float> gop(n, 10), gep(n, 0.2);
	
	adjust_gp(gop, gep, s);
	
	copy(gop.begin(), gop.end(), ostream_iterator<float>(cout, "\t"));
	cout << endl;
	exit(0);
}

int main(int argc, char* argv[])
{
	try
	{
		po::options_description desc("align options");
		desc.add_options()
			("help,h", "Display help message")
			("input,i", po::value<string>(), "Input file")
			("debug,d", po::value<int>(), "Debug output")
			("verbose,v", "Verbose output")
			("gap-open", po::value<float>(), "Gap open penalty")
			("gap-extend", po::value<float>(), "Gap extend penalty")
			("test,t", "run test function and exit")
			("guide-tree,g", po::value<string>(), "use existing guide tree")
			("matrix,m", po::value<string>(), "Substitution matrix, default is BLOSUM")
			;
	
		po::variables_map vm;
		po::store(po::parse_command_line(argc, argv, desc), vm);
		po::notify(vm);
		
		if (vm.count("help") or (vm.count("input") == 0 and vm.count("test") == 0))
		{
			cout << desc << endl;
			exit(1);
		}
		
		if (vm.count("debug"))
			DEBUG = vm["debug"].as<int>();
		VERBOSE = vm.count("verbose");

		if (vm.count("test"))
			test();
	
		// matrix
		string matrix = "PAM";
		if (vm.count("matrix"))
			matrix = vm["matrix"].as<string>();
		substitution_matrix_family mat(matrix);

		float gop = 10.f;	if (vm.count("gap-open")) gop = vm["gap-open"].as<float>();
		float gep = 0.2f;	if (vm.count("gap-extend")) gep = vm["gap-extend"].as<float>();

		fs::path path(vm["input"].as<string>());
		vector<entry> data;
		readFasta(path, data);
		
		if (data.size() < 2)
			throw my_bad("insufficient number of sequences");

		vector<entry*> alignment;
		base_node* root;
		
		if (data.size() == 2)
		{
			// no need to do difficult stuff, just align two sequences:
			
			float dist = calculateDistance(data[0], data[1]);
			root = new joined_node(new leaf_node(data[0]), new leaf_node(data[1]), dist / 2, dist / 2);
		}
		else
		{
			// create the leaf nodes
			vector<base_node*> tree;
			tree.reserve(data.size());
			foreach (entry& e, data)
				tree.push_back(new leaf_node(e));
			
			// calculate guide tree
			
			if (vm.count("guide-tree"))
				useGuideTree(vm["guide-tree"].as<string>(), tree);
			else
			{
				// a distance matrix
				distance_matrix<float> d(data.size());
				
				int n = 0;
				for (uint32 a = 0; a < data.size() - 1; ++a)
				{
					for (uint32 b = a + 1; b < data.size(); ++b)
					{
						if ((++n % 1000) == 0)
						{
							cerr << '.';
							if ((n % 60000) == 0)
								cerr << ' ' << n << endl;
						}
						
						d(a, b) = calculateDistance(data[a], data[b]);
					}
				}
				cerr << endl;

				joinNeighbours(d, tree);
			}
			
			root = tree.front();
		}
		
		if (VERBOSE)
			cout << *root << ';' << endl;
		
		createAlignment(root, alignment, mat, gop, gep);

		sort(alignment.begin(), alignment.end(),
			boost::bind(&entry::nr, _1) < boost::bind(&entry::nr, _2));

		report(alignment);

		delete root;
	}
	catch (const char* e)
	{
		cerr << e << endl;
		exit(1);
	}
	catch (const exception& e)
	{
		cerr << e.what() << endl;
		exit(1);
	}
	
	return 0;
}

