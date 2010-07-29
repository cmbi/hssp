// substitution matrix for multiple sequence alignments

#pragma once

#include "align.h"

#include <string>
#include <istream>
#include <cassert>

class substitution_matrix
{
  public:
						substitution_matrix(const std::string& name);
						substitution_matrix(
							const substitution_matrix& m, bool positive);

	virtual				~substitution_matrix() {}

	int32				operator()(aa a, aa b) const
						{
							return m_matrix(a, b);
						}

	float				mismatch_average() const		{ return m_mismatch_average; }

  private:
						substitution_matrix(
							const substitution_matrix&);
	substitution_matrix&
						operator=(const substitution_matrix&);

	void				read(std::istream& is);

	matrix<int8>		m_matrix;
	float				m_mismatch_average;
};

class substitution_matrix_family
{
  public:
						substitution_matrix_family(
							const std::string& name);

						~substitution_matrix_family();

	const substitution_matrix&
						operator()(float distance, bool positive) const
						{
							const substitution_matrix* result;
							
							if (positive)
							{
								if (distance >= 0.8f)
									result = m_pos_smat[0];
								else if (distance >= 0.6f)
									result = m_pos_smat[1];
								else if (distance >= 0.4f)
									result = m_pos_smat[2];
								else
									result = m_pos_smat[3];
							}
							else
							{
								if (distance >= 0.8f)
									result = m_smat[0];
								else if (distance >= 0.6f)
									result = m_smat[1];
								else if (distance >= 0.4f)
									result = m_smat[2];
								else
									result = m_smat[3];
							}
							
							return *result;
						}

  private:
						substitution_matrix_family(
							const substitution_matrix_family&);
	substitution_matrix_family&
						operator=(const substitution_matrix_family&);

	substitution_matrix*
						m_smat[4];
	substitution_matrix*
						m_pos_smat[4];
};

//ostream& operator<<(ostream& os, substitution_matrix& m)
//{
//	// print header
//	os << ' ';
//	for (uint32 i = 0; i < sizeof(kAA); ++i)
//		os << "  " << kAA[i];
//	os << endl;
//	
//	// print matrix
//	for (uint32 r = 0; r < sizeof(kAA); ++r)
//	{
//		os << kAA[r];
//		
//		for (uint32 c = 0; c < sizeof(kAA); ++c)
//			os << setw(3) << m(r, c);
//
//		os << endl;
//	}
//	
//	return os;
//}
