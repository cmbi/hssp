// substitution matrix code

#include "matrix.h"

#include <sstream>
#include <iostream>

#include <boost/iostreams/device/array.hpp>
#include <boost/iostreams/stream.hpp>
#include <boost/foreach.hpp>
#define foreach BOOST_FOREACH
#include <boost/bind.hpp>

#include "utils.h"
#include "../mtrx/matrices.h"

using namespace std;
namespace io = boost::iostreams;

substitution_matrix::substitution_matrix(const string& name)
	: m_matrix(sizeof(kAA), sizeof(kAA), 0)
	, m_scale_factor(0.75f)
{
	const MatrixInfo* mi = find_if(kMatrices, kMatrices + kMatrixCount,
		boost::bind(&MatrixInfo::name, _1) == name);
	
	if (mi == kMatrices + kMatrixCount)
		throw mas_exception(boost::format("missing matrix %1%") % name);
	
	io::stream<io::array_source> in(mi->m_data, strlen(mi->m_data));
	read(in);
}

substitution_matrix::substitution_matrix(
	const substitution_matrix& m, bool positive)
	: m_matrix(sizeof(kAA), sizeof(kAA), 0)
	, m_scale_factor(0.75f)
{
	int8 min = 0;
	
	for (uint32 y = 0; y < kAACount; ++y)
	{
		for (uint32 x = 0; x < kAACount; ++x)
		{
			m_matrix(x, y) = m.m_matrix(x, y);
			
			if (min > m_matrix(x, y))
				min = m_matrix(x, y);
		}
	}
	
	if (min < 0)
	{
		min = -min;
		
		for (uint32 y = 0; y < kAACount; ++y)
		{
			for (uint32 x = 0; x <= y; ++x)
			{
				m_matrix(x, y) += min;
				if (x != y)
					m_matrix(y, x) += min;
			}
		}
		
		float sum = 0;
		for (uint32 ry = 1; ry < 20; ++ry)
		{
			for (uint32 rx = 0; rx < ry; ++rx)
				sum += m_matrix(rx, ry);
		}
		
		m_mismatch_average = sum / ((20 * 19) / 2);
	}
}

void substitution_matrix::read(istream& is)
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
			throw mas_exception("invalid matrix file");
		
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
			if (v < numeric_limits<int8>::min() or v > numeric_limits<int8>::max())
				throw mas_exception("Invalid matrix, value out of range");

			m_matrix(row, ix[i]) = static_cast<int8>(v);
		}
	}
	
	// calculate mismatch_average
	
	float sum = 0;
	for (uint32 ry = 1; ry < 20; ++ry)
	{
		for (uint32 rx = 0; rx < ry; ++rx)
			sum += m_matrix(rx, ry);
	}
	
	m_mismatch_average = sum / ((20 * 19) / 2);
}

// --------------------------------------------------------------------

substitution_matrix_family::substitution_matrix_family(
	const std::string& name)
{
	if (name != "BLOSUM" and name != "PAM" and name != "GONNET")
		throw mas_exception(boost::format("unsuppported matrix %1%") % name);

	if (name == "BLOSUM")
	{
		m_cutoff[0] = 0.8f;
		m_smat[0] = new substitution_matrix(name + "80");
		m_cutoff[0] = 0.6f;
		m_smat[1] = new substitution_matrix(name + "62");
		m_cutoff[0] = 0.3f;
		m_smat[2] = new substitution_matrix(name + "45");
		m_cutoff[0] = 0;
		m_smat[3] = new substitution_matrix(name + "30");
	}
	else if (name == "PAM")
	{
		m_cutoff[0] = 0.8f;
		m_smat[0] = new substitution_matrix(name + "20");
		m_cutoff[0] = 0.6f;
		m_smat[1] = new substitution_matrix(name + "60");
		m_cutoff[0] = 0.4f;
		m_smat[2] = new substitution_matrix(name + "120");
		m_cutoff[0] = 0;
		m_smat[3] = new substitution_matrix(name + "350");
	}
	else //if (name == "GONNET")
	{
		m_cutoff[0] = 0.0f;
		m_smat[0] = new substitution_matrix(name + "250");
		m_smat[3] = m_smat[2] = m_smat[1] = NULL;
	}

	m_pos_smat[0] = new substitution_matrix(*m_smat[0], true);
	
	if (m_smat[1] != NULL)
		m_pos_smat[1] = new substitution_matrix(*m_smat[1], true);
	else
		m_pos_smat[1] = NULL;

	if (m_smat[2] != NULL)
		m_pos_smat[2] = new substitution_matrix(*m_smat[2], true);
	else
		m_pos_smat[2] = NULL;

	if (m_smat[3] != NULL)
		m_pos_smat[3] = new substitution_matrix(*m_smat[3], true);
	else
		m_pos_smat[3] = NULL;
}

substitution_matrix_family::~substitution_matrix_family()
{
	delete m_smat[0];
	delete m_smat[1];
	delete m_smat[2];
	delete m_smat[3];

	delete m_pos_smat[0];
	delete m_pos_smat[1];
	delete m_pos_smat[2];
	delete m_pos_smat[3];
}
