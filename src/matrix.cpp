// substitution matrix code

#include "matrix.h"

#include <boost/iostreams/device/array.hpp>
#include <boost/iostreams/stream.hpp>

#include "../matrices/matrices.h"

using namespace std;
namespace io = boost::iostreams;

class substitution_matrix_impl : public substitution_matrix
{
  public:
						substitution_matrix_impl(const string& name);

	void				read(istream& is);
};

substitution_matrix_impl::substitution_matrix_impl(const string& name)
{
	if (name == "BLOSUM80")
	{
		io::stream<io::array_source> in(kBLOSUM80, strlen(kBLOSUM80));
		read(in);
	}
	else if (name == "BLOSUM62")
	{
		io::stream<io::array_source> in(kBLOSUM62, strlen(kBLOSUM62));
		read(in);
	}
	else if (name == "BLOSUM45")
	{
		io::stream<io::array_source> in(kBLOSUM45, strlen(kBLOSUM45));
		read(in);
	}
	else if (name == "BLOSUM30")
	{
		io::stream<io::array_source> in(kBLOSUM30, strlen(kBLOSUM30));
		read(in);
	}
	else
		throw my_bad(boost::format("unsupported matrix %1%") % name);
}

void substitution_matrix_impl::read(istream& is)
{
	sequence ix;
	
	// first read up until we've got the header and calculate the index
	for (;;)
	{
		string line;
		getline(is, line);
		if (line.empty())
		{
			if (is.eof())
				break;
			continue;
		}
		if (line[0] == '#')
			continue;
		
		if (line[0] != ' ')
			throw my_bad("invalid matrix file");
		
		string h;
		foreach (char ch, line)
		{
			if (ch != ' ')
				h += ch;
		}
		
		ix = encode(h);
		
		break;
	}
	
	for (;;)
	{
		string line;
		getline(is, line);
		if (line.empty())
		{
			if (is.eof())
				break;
			continue;
		}
		if (line[0] == '#')
			continue;
		
		uint32 row = encode(line.substr(0, 1))[0];
		
		stringstream s(line.substr(1));
		int32 v;

		for (uint32 i = 0; i < ix.length(); ++i)
		{
			s >> v;
			m_matrix(row, ix[i]) = v;
		}
	}
}

// --------------------------------------------------------------------

substitution_matrix_family::substitution_matrix_family(
	const std::string& name)
{
	if (name != "BLOSUM")
		throw my_bad(boost::format("unsuppported matrix %1%") % name);

	m_smat[0] = new substitution_matrix_impl(name + "80");
	m_smat[1] = new substitution_matrix_impl(name + "62");
	m_smat[2] = new substitution_matrix_impl(name + "45");
	m_smat[3] = new substitution_matrix_impl(name + "30");
}

substitution_matrix_family::~substitution_matrix_family()
{
	delete m_smat[0];
	delete m_smat[1];
	delete m_smat[2];
	delete m_smat[3];
}
