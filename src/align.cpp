// align.cpp - simple attempt to write a multiple sequence alignment application
//

#include "align.h"

#include <iostream>
#include <iomanip>
#include <string>
#include <sstream>
#include <limits>

#include <boost/program_options.hpp>
#include <boost/format.hpp>
#include <boost/filesystem/fstream.hpp>
#include <boost/foreach.hpp>
#include <boost/tr1/tuple.hpp>
#define foreach BOOST_FOREACH
#include <boost/iostreams/device/array.hpp>
#include <boost/iostreams/stream.hpp>

#include "matrices.h"

using namespace std;
using namespace tr1;
namespace fs = boost::filesystem;
namespace po = boost::program_options;
namespace io = boost::iostreams;

int DEBUG = 0;

// --------------------------------------------------------------------

aa kAA_Reverse[256];

struct entry
{
					entry(const string& id, const sequence& seq)
						: m_id(id)
						, m_seq(seq)
						, m_weight(1.0f) {}

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
	
	uint32			dim() const			{ return m_n; }
	
//	tuple<uint32,uint32>
//					min() const;

	// erase two rows, add one at the end (for neighbour joining)
	void			erase_2(uint32 i, uint32 j);
	
	void			print(ostream& s);

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
T distance_matrix<T>::operator()(uint32 i, uint32 j) const
{
	assert(i < m_n); assert(j < m_n);
	
	T result = 0;

	if (i > j)
		result = m_data[(i * (i - 1)) / 2 + j];
	else if (i < j)
		result = m_data[(j * (j - 1)) / 2 + i];
	
	return result;
}

template<typename T>
inline
T& distance_matrix<T>::operator()(uint32 i, uint32 j)
{
	assert(i != j); assert(i < m_n); assert(j < m_n);
	
	if (i > j)
		return m_data[(i * (i - 1)) / 2 + j];
	else if (i < j)
		return m_data[(j * (j - 1)) / 2 + i];
}

//template<typename T>
//tuple<uint32,uint32> distance_matrix<T>::min() const
//{
//	uint32 min_i = 1, min_j = 0,
//		m = m_data[0], r = (m_n * (m_n - 1)) / 2;
//	
//	for (uint32 x = 1, i = 2, j = 0; x < r; ++x)
//	{
//		if (m_data[x] < m)
//		{
//			min_i = i;
//			min_j = j;
//			m = m_data[x];
//		}
//		
//		++j;
//		if (j == i)
//		{
//			j = 0;
//			++i;
//		}
//	}
//	
//	return make_tuple(min_i, min_j);
//}

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
ostream& operator<<(ostream& lhs, const distance_matrix<T>& rhs)
{
	for (uint32 i = 0; i < rhs.dim(); ++i)
	{
		for (uint32 j = 0; j < rhs.dim(); ++j)
		{
			if (i == j)
				lhs << '-';
			else
				lhs << rhs(i, j);
			
			if (j < rhs.dim() - 1)
				lhs << '\t';
		}
		lhs << endl;
	}

	return lhs;
}

struct trace_back
{
				trace_back(
					uint32		inDimX,
					uint32		inDimY)
					: mDimX(inDimX)
					, mDimY(inDimY)
					, mData(NULL)
				{
					mData = new int8[(mDimX + 1) * (mDimY + 1)];
					memset(mData, 0, (mDimX + 1) * (mDimY + 1));
				}
				
				~trace_back()
				{
					delete[] mData;
				}
	
	int32		operator()(
					int32		inB,
					int32		inIx,
					int32		inIy,
					uint32		inX,
					uint32		inY)
				{
					int32 result;

					if (inB >= inIx and inB >= inIy)
					{
						result = inB;
						set(inX, inY, 0);
					}
					else if (inIx >= inB and inIx >= inIy)
					{
						result = inIx;
						set(inX, inY, 1);
					}
					else
					{
						result = inIy;
						set(inX, inY, -1);
					}
					
					return result;
				}

	void		set(
					uint32		inX,
					uint32		inY,
					int16		inV)
				{
					assert(inX <= mDimX);
					assert(inY <= mDimY);
					mData[inX * mDimY + inY] = inV;
				}
	
	int16		test(
					uint32		inX,
					uint32		inY)
				{
					assert(inX <= mDimX);
					assert(inY <= mDimY);
					return mData[inX * mDimY + inY];
				}

	int16		operator()(
					uint32		inX,
					uint32		inY) const
				{
					assert(inX <= mDimX);
					assert(inY <= mDimY);
					return mData[inX * mDimY + inY];
				}

	uint32		mDimX, mDimY;
	int8*		mData;
};

// --------------------------------------------------------------------



// --------------------------------------------------------------------

class pss_matrix
{
  public:
						pss_matrix(const substitution_matrix& mat,
							const vector<entry>& alignment,
							int16 gapOpen, int16 gapExtend);

	int32				operator()(uint32 p, aa a) const
						{
							return m_matrix[p].m_scores[a];
						}
						
	int32				gap_open(uint32 p) const
						{
							return m_matrix[p].m_gap_open;
						}

	int32				gap_extend(uint32 p) const
						{
							return m_matrix[p].m_gap_extend;
						}

  private:

	struct ps_score {
		int16			m_scores[kAA_Count];
		int16			m_gap_open;
		int16			m_gap_extend;
	};

	vector<ps_score>	m_matrix;
};

pss_matrix::pss_matrix(const substitution_matrix& mat,
	const vector<entry>& alignment, int16 gapOpen, int16 gapExtend)
{
	uint32 nseq = alignment.size();
	assert(nseq > 0);
	
	uint32 slen = alignment[0].m_seq.length();
	m_matrix.reserve(slen);
	
	for (uint32 p = 0; p < slen; ++p)
	{
		int32 s[kAA_Count] = { };
		int32 n = 0, gaps = 0;

		for (uint32 seq = 0; seq < nseq; ++seq)
		{
			aa r = alignment[seq].m_seq[p];
			if (r == kSignalGapCode)
				++gaps;
			else if (r < kFilteredCode)
			{
				++n;
				for (uint32 a = 0; a < kAA_Count; ++a)
					s[a] += mat(a, r);
			}
		}
		
		if (n > 1)
		{
			for (uint32 a = 0; a < kAA_Count; ++a)
				s[a] /= n;
		}
	
		ps_score ms;
		for (uint32 aa = 0; aa < kAA_Count; ++aa)
		{
			if (s[aa] < kSentinalScore)
				ms.m_scores[aa] = kSentinalScore;
			else
				ms.m_scores[aa] = static_cast<int16>(s[aa]);
		}
		
		ms.m_gap_open = gapOpen;
		ms.m_gap_extend = gapExtend;

//		if (gaps > 0)
//		{
//			assert(gaps < static_cast<int32>(nseq));
//			
//			gaps *= inGapScaleFactor;
//			
//			float gapfactor = static_cast<float>(nseq - gaps) / nseq;
//			
//			ms.gap_open *= gapfactor;
//			ms.gap_extend *= gapfactor;
//		}
		
		m_matrix.push_back(ms);
	}
}

// --------------------------------------------------------------------

struct base_node
{
	virtual				~base_node() {}

	virtual void		print(ostream& s) = 0;	

	virtual base_node*	left() const		{ return 0; }
	virtual base_node*	right() const		{ return 0; }
};

struct joined_node : public base_node
{
						joined_node(base_node* left, base_node* right,
							int32 d_left, int32 d_right)
							: m_left(left)
							, m_right(right)
							, m_d_left(d_left)
							, m_d_right(d_right) {}

	virtual				~joined_node()
						{
							delete m_left;
							delete m_right;
						}

	virtual void		print(ostream& s)
						{
							s << '(' << endl;
							m_left->print(s);
							s << ':' << m_d_left << ',' << endl;
							m_right->print(s);
							s << ':' << m_d_right << endl
							  << ')';
						}

	virtual base_node*	left() const		{ return m_left; }
	virtual base_node*	right() const		{ return m_right; }
	
	base_node*			m_left;
	base_node*			m_right;
	float				m_d_left;
	float				m_d_right;
	uint32				m_leaf_count;
};

struct leaf_node : public base_node
{
						leaf_node(const entry& e)
							: m_name(e.m_id)
							, m_sequence(e.m_seq) {}

	virtual void		print(ostream& s)
						{
							s << m_name;
						}

	string				m_name;
	sequence			m_sequence;
};

ostream& operator<<(ostream& lhs, base_node* rhs)
{
	rhs->print(lhs);
	lhs << ';' << endl;
	return lhs;
}

// --------------------------------------------------------------------
// compute the distance between two sequences using the
// Levenshtein algorithm.

uint16 calculateDistance(const sequence& s, const sequence& t)
{
	uint32 m = s.length();
	uint32 n = t.length();
	
	matrix<uint16> d(m + 1, n + 1);
	
	for (uint16 i = 0; i <= m; ++i)
		d(i, 0) = i;
	for (uint16 j = 0; j <= n; ++j)
		d(0, j) = j;
	
	for (uint16 j = 1; j <= n; ++j)
	{
		for (uint16 i = 1; i <= m; ++i)
		{
			if (s[i - 1] == t[j - 1])
				d(i, j) = d(i - 1, j - 1);
			else
			{
				uint16 v1 = d(i - 1, j) + 1;
				uint16 v2 = d(i, j - 1) + 1;
				uint16 v3 = d(i - 1, j - 1) + 1;
				
				if (v1 > v2)
					v1 = v2;
				if (v1 > v3)
					v1 = v3;
				
				d(i, j) = v1;
			}
		}
	}
	
	return d(m, n);
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
				seq.push_back(entry(id, encode(s)));
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

void joinNeighbours(distance_matrix<uint16>& d, vector<base_node*>& tree)
{
	uint32 r = tree.size();
	
//	matrix<int32> Q(r);
	vector<int32> sum(r);
	
	// calculate the sums first
	for (uint32 i = 0; i < r; ++i)
	{
		for (uint32 j = 0; j < r; ++j)
		{
			if (i != j)
				sum[i] += d(i, j);
		}
	}
	
	// calculate Q
	uint32 min_i = 0, min_j = 0;
	int32 m = numeric_limits<int32>::max();
	
	for (uint32 i = 1; i < r; ++i)
	{
		for (uint32 j = 0; j < i; ++j)
		{
			int32 v = d(i, j) - abs(sum[i] + sum[j]) / (r - 2);

			if (m > v)
			{
				min_i = i;
				min_j = j;
				m = v;
			}
//			
//			Q(i, j) = v;
		}
	}
	
//	// debug output
//	cerr << "Q is now:" << endl << Q << endl
//		 << "joining " << min_i << " and " << min_j << endl;

	// distance to joined node
	int32 d_i, d_j;
	int32 half_dij = d(min_i, min_j) / 2;
	
	d_i = half_dij + abs(sum[min_i] - sum[min_j]) / (2 * (r - 2));
	d_j = d(min_i, min_j) - d_i;

	joined_node* jn = new joined_node(tree[min_i], tree[min_j], d_i, d_j);
	assert(min_j < min_i);
	tree.erase(tree.begin() + min_i);
	tree.erase(tree.begin() + min_j);
	tree.push_back(jn);
	
	// distances to other nodes
	vector<int32> dn; dn.reserve(r - 2);
	for (uint32 x = 0; x < r; ++x)
	{
		if (x == min_i or x == min_j)
			continue;
		dn.push_back(d(x, min_i) + d(x, min_j) - half_dij);
	}
	
	// fill new distance matrix
	d.erase_2(min_i, min_j);
	--r;
	for (uint32 x = 0; x < r - 1; ++x)
		d(x, r - 1) = dn[x];
}

// --------------------------------------------------------------------

int32 score(const vector<entry>& a, const vector<entry>& b,
	uint32 ix_a, uint32 ix_b, const substitution_matrix& mat)
{
	float result = 0;
	
	foreach (const entry& ea, a)
	{
		foreach (const entry &eb, b)
		{
			assert(ix_a < ea.m_seq.length());
			assert(ix_b < eb.m_seq.length());

			result += ea.m_weight * eb.m_weight * mat(ea.m_seq[ix_a], eb.m_seq[ix_b]);
		}
	}
	
	return static_cast<int32>(result / (a.size() * b.size()));
}

void align(vector<entry>& a, vector<entry>& b, vector<entry>& c,
	const substitution_matrix& mat,
	int32 gapOpen = 10, int32 gapExtend = 1)
{
	int32 dimX = a.front().m_seq.length();
	int32 dimY = b.front().m_seq.length();
	
	matrix<int32> B(dimX + 1, dimY + 1);
	matrix<int32> Ix(dimX + 1, dimY + 1);
	matrix<int32> Iy(dimX + 1, dimY + 1);
	
	Ix(0, 0) = 0;
	Ix(0, 1) = 0;
	Iy(0, 0) = 0;
	Iy(1, 0) = 0;
	B(0, 0) = 0;

	trace_back tb(dimX, dimY);

	for (int32 x = 1; x <= dimX; ++x)
	{
		for (int32 y = 1; y <= dimY; ++y)
		{
			int32 Ix1 = Ix(x - 1, y);	if (x == 1 and y > 1) Ix1 = -(gapOpen + y * gapExtend);
			int32 Iy1 = Iy(x, y - 1);	if (x > 1 and y == 1) Iy1 = -(gapOpen + x * gapExtend);
			
			// (1)
			int32 M = score(a, b, x - 1, y - 1, mat);
			if (x > 1 and y > 1)
				M += B(x - 1, y - 1);

			if (M >= Ix1 and M >= Iy1)
			{
				tb.set(x, y, 0);
				B(x, y) = M;
			}
			else if (Ix1 >= Iy1)
			{
				tb.set(x, y, 1);
				B(x, y) = Ix1;
			}
			else
			{
				tb.set(x, y, -1);
				B(x, y) = Iy1;
			}

if (DEBUG)
{
cout << "x: " << x << "; y: " << y << "; m: " << M << "; Ix1: " << Ix1 << "; Iy1: " << Iy1
	 << "; B=> " << B(x, y) << "; tb=> " << tb(x, y) << endl;
}
//			int32 d = s.GapOpen(y - 1);
//			int32 e = s.GapExtend(y - 1);
			int32 d = gapOpen;
			int32 e = gapExtend;

			// (3)
			Ix(x, y) = max(M - d, Ix1 - e);
			
			// (4)
			Iy(x, y) = max(M - d, Iy1 - e);
		}
	}
	
	// build the final alignment
	
	int32 x = dimX, y = dimY;
	while (x > 0 and y > 0)
	{
		if (tb(x, y) == 0)
		{
			--x;
			--y;
		}
		else if (tb(x, y) < 0)
		{
			foreach (entry& e, a)
				e.m_seq.insert(e.m_seq.begin() + x, char(kSignalGapCode));
			
			--y;
		}
		else
		{
			foreach (entry& e, b)
				e.m_seq.insert(e.m_seq.begin() + y, char(kSignalGapCode));
			
			--x;
		}
	}

	while (x > 0)
	{
		foreach (entry& e, b)
			e.m_seq.insert(e.m_seq.begin() + y, char(kSignalGapCode));
		
		--x;
	}
	
	while (y > 0)
	{
		foreach (entry& e, a)
			e.m_seq.insert(e.m_seq.begin() + x, char(kSignalGapCode));
		
		--y;
	}
			

	copy(a.begin(), a.end(), back_inserter(c));
	copy(b.begin(), b.end(), back_inserter(c));
	
// debug code:
	if (DEBUG)
	{
		foreach (entry& e, c)
			cout << e.m_id << ": " << decode(e.m_seq) << endl;
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

void alignAlignments(vector<entry>& a, vector<entry>& b, vector<entry>& c,
	const substitution_matrix& mat)
{
	align(a, b, c, mat);
}

void createAlignment(base_node* node, vector<entry>& alignment,
	const substitution_matrix& mat)
{
	if (node->left() and node->right())
	{
		vector<entry> a, b;	
		createAlignment(node->left(), a, mat);
		createAlignment(node->right(), b, mat);
		
		alignAlignments(a, b, alignment, mat);
	}
	else
	{
		leaf_node* n = dynamic_cast<leaf_node*>(node);
		if (n == 0)
			throw my_bad("internal error (not a leaf node)");
		
		alignment.push_back(entry(n->m_name, n->m_sequence));
	}
}

// --------------------------------------------------------------------

void report(const vector<entry>& alignment, const substitution_matrix& mat)
{
	cout << "CLUSTAL FORMAT for MaartensAlignment" << endl;

	uint32 nseq = alignment.size();
	uint32 len = alignment[0].m_seq.length();
	uint32 offset = 0;
	while (offset < len)
	{
		uint32 n = alignment[0].m_seq.length() - offset;
		if (n > 60)
			n = 60;
		
		vector<set<aa> > dist(n);
		
		foreach (const entry& e, alignment)
		{
			sequence ss = e.m_seq.substr(offset, n);
			
			for (uint32 i = 0; i < n; ++i)
				dist[i].insert(ss[i]);

			string id = e.m_id;
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

int main(int argc, char* argv[])
{
	try
	{
		po::options_description desc("align options");
		desc.add_options()
			("help,h", "Display help message")
			("input,i", po::value<string>(), "Input file")
			("debug,d", "Debug output")
			("matrix,m", po::value<string>(), "Substitution matrix, default is BLOSUM")
			;
	
		po::variables_map vm;
		po::store(po::parse_command_line(argc, argv, desc), vm);
		po::notify(vm);
		
		if (vm.count("help") or vm.count("input") == 0)
		{
			cout << desc << endl;
			exit(1);
		}
	
		DEBUG = vm.count("debug");
	
		// matrix
		string matrix = "BLOSUM";
		if (vm.count("matrix"))
			matrix = vm["matrix"].as<string>();
		substitution_matrix mat(matrix);
//		cout << "matrix:" << endl << mat << endl;

		fs::path path(vm["input"].as<string>());
		vector<entry> data;
		readFasta(path, data);
		
		if (data.size() < 2)
			throw my_bad("insufficient number of sequences");
		
		// a distance matrix
		distance_matrix<uint16> d(data.size());
		
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
				
				d(a, b) = calculateDistance(data[a].m_seq, data[b].m_seq);
			}
		}
		cerr << endl;
		
		// create the leaf nodes
		vector<base_node*> tree;
		tree.reserve(data.size());
		foreach (const entry& e, data)
			tree.push_back(new leaf_node(e));
		
		// calculate guide tree
		while (tree.size() > 2)
			joinNeighbours(d, tree);
		
		joined_node* root = new joined_node(tree[0], tree[1], d(0, 1) / 2, d(0, 1) / 2);
		
//		cout << root << endl;
		
		vector<entry> alignment;
		createAlignment(root, alignment, mat);
		report(alignment, mat);

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

