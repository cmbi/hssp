// align.cpp - simple attempt to write a multiple sequence alignment application
//

typedef short			int16;
typedef unsigned short	uint16;
typedef long			int32;
typedef unsigned long	uint32;

#include <iostream>
#include <string>
#include <limits>

#include <boost/program_options.hpp>
#include <boost/format.hpp>
#include <boost/filesystem/fstream.hpp>
#include <boost/foreach.hpp>
#include <boost/tr1/tuple.hpp>
#define foreach BOOST_FOREACH

using namespace std;
using namespace tr1;
namespace fs = boost::filesystem;
namespace po = boost::program_options;

#if defined(_MSC_VER)
#include <ciso646>
#define snprintf _snprintf
#endif

class my_bad : public exception
{
  public:
					my_bad(const char* msg)
					{
						snprintf(m_msg, sizeof(m_msg), "%s", msg);
					}

					my_bad(const string& msg)
					{
						snprintf(m_msg, sizeof(m_msg), "%s", msg.c_str());
					}

					my_bad(const boost::format& msg)
					{
						snprintf(m_msg, sizeof(m_msg), "%s", msg.str().c_str());
					}

	virtual const char*
					what() const throw()	{ return m_msg; }

  private:
	char			m_msg[1024];
};

template<typename T>
class matrix
{
  public:
	typedef T value_type;
	
					matrix(uint32 m, uint32 n);
	virtual			~matrix();
	
	value_type		operator()(uint32 i, uint32 j) const;
	value_type&		operator()(uint32 i, uint32 j);

  private:
	value_type*		m_data;
	uint32			m_m, m_n;
};

template<typename T>
matrix<T>::matrix(uint32 m, uint32 n)
	: m_m(m)
	, m_n(n)
{
	m_data = new value_type[m_m * m_n];
}

template<typename T>
matrix<T>::~matrix()
{
	delete[] m_data;
}

template<typename T>
inline
T matrix<T>::operator()(uint32 i, uint32 j) const
{
	assert(i < m_m); assert(j < m_n);
	return m_data[i + j * m_m];
}

template<typename T>
inline
T& matrix<T>::operator()(uint32 i, uint32 j)
{
	assert(i < m_m); assert(j < m_n);
	return m_data[i + j * m_m];
}

template<typename T>
class qmatrix
{
  public:
	typedef T value_type;
	
					qmatrix(uint32 n);
	virtual			~qmatrix();
	
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
qmatrix<T>::qmatrix(uint32 n)
	: m_n(n)
{
	m_data = new value_type[(m_n * (m_n - 1)) / 2];
}

template<typename T>
qmatrix<T>::~qmatrix()
{
	delete[] m_data;
}

template<typename T>
inline
T qmatrix<T>::operator()(uint32 i, uint32 j) const
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
T& qmatrix<T>::operator()(uint32 i, uint32 j)
{
	assert(i != j); assert(i < m_n); assert(j < m_n);
	
	if (i > j)
		return m_data[(i * (i - 1)) / 2 + j];
	else if (i < j)
		return m_data[(j * (j - 1)) / 2 + i];
}

//template<typename T>
//tuple<uint32,uint32> qmatrix<T>::min() const
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
void qmatrix<T>::erase_2(uint32 di, uint32 dj)
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
ostream& operator<<(ostream& lhs, const qmatrix<T>& rhs)
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
	int32				m_d_left;
	int32				m_d_right;
};

struct leaf_node : public base_node
{
						leaf_node(const string& name)
							: m_name(name) {}

	virtual void		print(ostream& s)
						{
							s << m_name;
						}

	string				m_name;
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

uint16 calculateDistance(const string& s, const string& t)
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
void readFasta(fs::path path, vector<pair<string,string> >& seq)
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
				seq.push_back(make_pair(id, s));
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

void joinNeighbours(qmatrix<uint16>& d, vector<base_node*>& tree)
{
	uint32 r = tree.size();
	
//	qmatrix<int32> Q(r);
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

int main(int argc, char* argv[])
{
	try
	{
		po::options_description desc("align options");
		desc.add_options()
			("help,h", "Display help message")
			("input,i", po::value<string>(), "Input file")
			;
	
		po::variables_map vm;
		po::store(po::parse_command_line(argc, argv, desc), vm);
		po::notify(vm);
		
		if (vm.count("help") or vm.count("input") == 0)
		{
			cout << desc << endl;
			exit(1);
		}
	
		fs::path path(vm["input"].as<string>());
		
		typedef pair<string,string> entry;
		vector<entry> seq;
		readFasta(path, seq);
		
		if (seq.size() < 2)
			throw my_bad("insufficient number of sequences");
		
		// a distance matrix
		qmatrix<uint16> d(seq.size());
		
		int n = 0;
		for (uint32 a = 0; a < seq.size() - 1; ++a)
		{
			for (uint32 b = a + 1; b < seq.size(); ++b)
			{
				if ((++n % 1000) == 0)
				{
					cerr << '.';
					if ((n % 60000) == 0)
						cerr << ' ' << n << endl;
				}
				
				d(a, b) = calculateDistance(seq[a].second, seq[b].second);
			}
		}
		cerr << endl;
		
		// build 'tree'
		vector<base_node*> tree;
		tree.reserve(seq.size());
		foreach (const entry& e, seq)
			tree.push_back(new leaf_node(e.first));
		
		// calculate initial Q
		while (tree.size() > 2)
			joinNeighbours(d, tree);
		
		joined_node* root = new joined_node(tree[0], tree[1], d(0, 1) / 2, d(0, 1) / 2);
		
		cout << root << endl;

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

