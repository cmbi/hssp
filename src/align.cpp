// align.cpp - simple attempt to write a multiple sequence alignment application
//

typedef short			int16;
typedef unsigned short	uint16;
typedef long			int32;
typedef unsigned long	uint32;

#include <iostream>
#include <string>

#include <boost/program_options.hpp>
#include <boost/format.hpp>
#include <boost/filesystem/fstream.hpp>
#include <boost/foreach.hpp>
#define foreach BOOST_FOREACH

using namespace std;
namespace fs = boost::filesystem;
namespace po = boost::program_options;

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
	return m_data[i + j * m_m];
}

template<typename T>
inline
T& matrix<T>::operator()(uint32 i, uint32 j)
{
	return m_data[i + j * m_m];
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
		throw my_bad((boost::format("input file '%1%' not opened") % path.string()).str());
	
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
		
//		foreach (const entry& e, seq)
//			cout << '>' << e.first << endl << e.second << endl;
		
		for (uint32 a = 0; a < seq.size() - 1; ++a)
		{
			for (uint32 b = a + 1; b < seq.size(); ++b)
			{
				uint16 d = calculateDistance(seq[a].second, seq[b].second);
				cout << seq[a].first << '\t' << seq[b].first << '\t' << d << endl;
			}
		}
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

