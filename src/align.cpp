// align.cpp - simple attempt to write a multiple sequence alignment application
//

#include "align.h"

#include <iostream>
#include <iomanip>
#include <string>
#include <sstream>
#include <limits>
#include <cmath>

#include <boost/program_options.hpp>
#include <boost/format.hpp>
#include <boost/filesystem/fstream.hpp>
#include <boost/foreach.hpp>
#include <boost/tr1/tuple.hpp>
#define foreach BOOST_FOREACH
#include <boost/iostreams/device/array.hpp>
#include <boost/iostreams/stream.hpp>
#include <boost/bind.hpp>

#include "matrix.h"

using namespace std;
using namespace tr1;
namespace fs = boost::filesystem;
namespace po = boost::program_options;
namespace io = boost::iostreams;

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

	uint32			m_nr;
	string			m_id;
	sequence		m_seq;
	float			m_weight;
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
	uint32				m_left_leaf_count, m_right_leaf_count, m_leaf_count;
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

joined_node::joined_node(base_node* left, base_node* right,
	float d_left, float d_right)
	: m_left(left)
	, m_right(right)
	, m_d_left(d_left)
	, m_d_right(d_right)
	, m_left_leaf_count(left->leaf_count())
	, m_right_leaf_count(right->leaf_count())
	, m_leaf_count(m_left_leaf_count + m_right_leaf_count)
{
	left->add_weight(d_left);
	right->add_weight(d_right);
}

void joined_node::add_weight(float w)
{
	w *= m_left_leaf_count;
	m_left->add_weight(w / m_left_leaf_count);
	m_right->add_weight(w / m_right_leaf_count);
}

// --------------------------------------------------------------------
// compute the distance between two sequences using the
// Levenshtein algorithm.

//uint16 calculateDistance(const sequence& s, const sequence& t)
//{
//	uint32 m = s.length();
//	uint32 n = t.length();
//	
//	matrix<uint16> d(m + 1, n + 1);
//	
//	for (uint16 i = 0; i <= m; ++i)
//		d(i, 0) = i;
//	for (uint16 j = 0; j <= n; ++j)
//		d(0, j) = j;
//	
//	for (uint16 j = 1; j <= n; ++j)
//	{
//		for (uint16 i = 1; i <= m; ++i)
//		{
//			if (s[i - 1] == t[j - 1])
//				d(i, j) = d(i - 1, j - 1);
//			else
//			{
//				uint16 v1 = d(i - 1, j) + 1;
//				uint16 v2 = d(i, j - 1) + 1;
//				uint16 v3 = d(i - 1, j - 1) + 1;
//				
//				if (v1 > v2)
//					v1 = v2;
//				if (v1 > v3)
//					v1 = v3;
//				
//				d(i, j) = v1;
//			}
//		}
//	}
//	
//	return d(m, n);
//}

float calculateDistance(const entry& a, const entry& b)
{
	int32 x, dimX = a.m_seq.length();
	int32 y, dimY = b.m_seq.length();
	
	matrix<float>	B(dimX + 1, dimY + 1);
	matrix<float>	Ix(dimX + 1, dimY + 1);
	matrix<float>	Iy(dimX + 1, dimY + 1);
	matrix<int8>	tb(dimX + 1, dimY + 1);
	
	Ix(0, 1) = 0;
	Iy(1, 0) = 0;

	static const substitution_matrix smat("GONNET250");
	static const float gapOpen = 10;
	static const float gapExtend = 0.1;

	for (x = 1; x <= dimX; ++x)
	{
		for (y = 1; y <= dimY; ++y)
		{
//			float Ix1 = Ix(x - 1, y);	if (x == 1 and y > 1) Ix1 = -(gapOpen + y * gapExtend);
//			float Iy1 = Iy(x, y - 1);	if (x > 1 and y == 1) Iy1 = -(gapOpen + x * gapExtend);
			float Ix1 = Ix(x - 1, y);	if (x == 1 and y > 1) Ix1 = 0;
			float Iy1 = Iy(x, y - 1);	if (x > 1 and y == 1) Iy1 = 0;
			
			// (1)
			float M = smat(a.m_seq[x], b.m_seq[y]);
			if (x > 1 and y > 1)
				M += B(x - 1, y - 1);

			if (M >= Ix1 and M >= Iy1)
			{
				tb(x, y) = 0;
				B(x, y) = M;
			}
			else if (Ix1 >= Iy1)
			{
				tb(x, y) = 1;
				B(x, y) = Ix1;
			}
			else
			{
				tb(x, y) = -1;
				B(x, y) = Iy1;
			}

			// (3)
			Ix(x, y) = max(M - gapOpen, Ix1 - gapExtend);
			
			// (4)
			Iy(x, y) = max(M - gapOpen, Iy1 - gapExtend);
		}
	}
	
	// traceback, counting the identities
	uint32 n = 0;
	x = dimX; y = dimY;
	while (x > 0 or y > 0)
	{
		if (x == 0 or tb(x, y) < 0)
			--y;
		else if (y == 0 or tb(x, y) > 0)
			--x;
		else
		{
			if (a.m_seq[x] == b.m_seq[y])
				++n;
			--x;
			--y;
		}
	}
	
	float result = 1.0f - float(n) / max(dimX, dimY);
	
	if (VERBOSE)
	{
		cout << "Sequences (" << a.m_nr + 1 << ':' << b.m_nr + 1 << ") Aligned. Score: " << setw(4) << setprecision(2) << result
//			 << "  score " << int32(B(dimX, dimY))
			 << endl;
	}
	
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

int32 score(const vector<entry*>& a, const vector<entry*>& b,
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
			
			float s;
			
			if (ra != kSignalGapCode and rb != kSignalGapCode)
				s = mat(ra, rb);

			result += ea->m_weight * eb->m_weight * s;
		}
	}
	
	return static_cast<int32>(result * 1.0f / (a.size() * b.size()));
}

void adjust_gp(vector<float>& gop, vector<float>& gep, const vector<entry*>& seq)
{
	assert(gop.size() == seq.front()->m_seq.length());

	vector<uint32> gaps(gop.size());

	foreach (const entry* e, seq)
	{
		for (uint32 ix = 0; ix < gop.size(); ++ix)
		{
			if (e->m_seq[ix] == kSignalGapCode)
				gaps[ix] += 1;
		}
	}
	
	for (int32 ix = 0; ix < static_cast<int32>(gop.size()); ++ix)
	{
		// if there is a gap, lower gap open cost
		if (gaps[ix] > 0)
		{
			gop[ix] *= 0.3f * (seq.size() - gaps[ix]);
			gep[ix] /= 2;
		}
		
		// else if there is a gap within 8 residues, increase gap cost
		else
		{
			for (int32 d = 0; d < 8; ++d)
			{
				if ((ix + d < gaps.size() and gaps[ix + d] > 0) or
					(ix - d >= 0 and gaps[ix - d] > 0))
				{
					cout << "ix: " << ix << " d: " << d << " f: " << (2 + abs((8 - d) * 2)) / 8.f << endl;
					
					gop[ix] += gop[ix] * (2 + abs((8 - d) * 2)) / 8.f;
					break;
				}
			}
		}
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
	
	matrix<float> B(dimX + 1, dimY + 1);
	matrix<float> Ix(dimX + 1, dimY + 1);
	matrix<float> Iy(dimX + 1, dimY + 1);
	matrix<int8> tb(dimX + 1, dimY + 1);
	
	Ix(0, 1) = 0;
	Iy(1, 0) = 0;

	const substitution_matrix& smat = mat_fam(abs(node->m_d_left + node->m_d_right), true);
	
	// initial gap costs
	gop = abs((gop + log(min(dimX, dimY))) * smat.mismatch_average()) / 10;
	gep = gep * (1 + abs(log(float(dimX) / dimY)));
	
	// position specific gap open costs
	vector<float> gop_a(dimX, gop), gep_a(dimX, gep);
	adjust_gp(gop_a, gep_a, a);
	
	vector<float> gop_b(dimY, gop), gep_b(dimY, gep);
	adjust_gp(gop_b, gep_b, b);

	if (VERBOSE)
	{
		cout << "gop: " << gop << "; "
			 << "gep: " << gep << endl;
		
		copy(gop_a.begin(), gop_a.end(), ostream_iterator<float>(cout, ", ")); cout << endl;
		copy(gop_b.begin(), gop_b.end(), ostream_iterator<float>(cout, ", ")); cout << endl;
	}

	for (x = 1; x <= dimX; ++x)
	{
		for (y = 1; y <= dimY; ++y)
		{
			float Ix1 = Ix(x - 1, y);	if (x == 1 and y > 1) Ix1 = -(gop + y * gep);
			float Iy1 = Iy(x, y - 1);	if (x > 1 and y == 1) Iy1 = -(gop + x * gep);
			
			// (1)
			float M = score(a, b, x - 1, y - 1, smat);
			if (x > 1 and y > 1)
				M += B(x - 1, y - 1);

			if (M >= Ix1 and M >= Iy1)
			{
				tb(x, y) = 0;
				B(x, y) = M;
			}
			else if (Ix1 >= Iy1)
			{
				tb(x, y) = 1;
				B(x, y) = Ix1;
			}
			else
			{
				tb(x, y) = -1;
				B(x, y) = Iy1;
			}

if (DEBUG)
{
	cout << "x: " << x << "; y: " << y << "; m: " << M << "; Ix1: " << Ix1 << "; Iy1: " << Iy1
		 << "; B=> " << B(x, y) << "; tb=> " << tb(x, y) << endl;
}

			// (3)
			Ix(x, y) = max(M - gop_a[x], Ix1 - gep_a[x]);
			
			// (4)
			Iy(x, y) = max(M - gop_b[y], Iy1 - gep_b[y]);
		}
	}
	
	// build the final alignment
	
	x = dimX;
	y = dimY;

	while (x > 0 or y > 0)
	{
		if (x == 0 or tb(x, y) < 0)
		{
			foreach (entry* e, a)
				e->m_seq.insert(e->m_seq.begin() + x, char(kSignalGapCode));
			--y;
		}
		else if (y == 0 or tb(x, y) > 0)
		{
			foreach (entry* e, b)
				e->m_seq.insert(e->m_seq.begin() + y, char(kSignalGapCode));
			--x;
		}
		else
		{
			assert(x > 0); assert(y > 0);
			--x;
			--y;
		}
	}

	copy(a.begin(), a.end(), back_inserter(c));
	copy(b.begin(), b.end(), back_inserter(c));
	
// debug code:
	if (DEBUG)
	{
		foreach (entry* e, c)
			cout << e->m_id << ": " << decode(e->m_seq) << endl;
		cout << endl;
		
		cout << "      ";
		for (int32 x = 1; x <= dimX; ++x)
			cout << setw(8) << x << "   ";
		cout << endl;
		
		for (int32 y = 1; y <= dimY; ++y)
		{
			cout << setw(6) << y;
			for (int32 x = 1; x <= dimX; ++x)
				cout << ' ' << setw(7) << B(x, y) << ' ' << setw(2) << tb(x, y);
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
		
		align(static_cast<joined_node*>(node), a, b, alignment, mat, gop, gep);
	}
	else
	{
		leaf_node* n = dynamic_cast<leaf_node*>(node);
		if (n == 0)
			throw my_bad("internal error (not a leaf node)");
		
		alignment.push_back(&n->m_entry);
	}
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
		
		vector<set<aa> > dist(n);
		
		foreach (const entry* e, alignment)
		{
			sequence ss = e->m_seq.substr(offset, n);
			
			for (uint32 i = 0; i < n; ++i)
				dist[i].insert(ss[i]);

			string id = e->m_id;
			if (id.length() > 15)
				id = id.substr(0, 12) + "...";
			else if (id.length() < 15)
				id += string(15 - id.length(), ' ');
			
			cout << id << ' ' << decode(ss) << endl;
		}
		
//		string scores(n, '*');
//		for (uint32 i = 0; i < n; ++i)
//		{
//			if (dist[i].size() > 1)
//			{
//				
//			}
//		}
//		
//		cout << string(16, ' ') << scores << endl;
		
		offset += n;
		cout << endl;
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
			("debug,d", "Debug output")
			("verbose,v", "Verbose output")
			("gap_open", po::value<float>(), "Gap open penalty")
			("gap_extend", po::value<float>(), "Gap extend penalty")
			("test,t", "run test function and exit")
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
		
		DEBUG = vm.count("debug");
		VERBOSE = vm.count("verbose");

		if (vm.count("test"))
			test();
	
		// matrix
		string matrix = "BLOSUM";
		if (vm.count("matrix"))
			matrix = vm["matrix"].as<string>();
		substitution_matrix_family mat(matrix);

		float gop = 10.f;	if (vm.count("gap_open")) gop = vm["gap_open"].as<float>();
		float gep = 0.2f;	if (vm.count("gap_extend")) gep = vm["gap_extend"].as<float>();

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
			
			// create the leaf nodes
			vector<base_node*> tree;
			tree.reserve(data.size());
			foreach (entry& e, data)
				tree.push_back(new leaf_node(e));
			
			// calculate guide tree
			
			joinNeighbours(d, tree);
			root = tree.front();
		}
		
		if (VERBOSE)
			cout << *root << endl;
		
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

