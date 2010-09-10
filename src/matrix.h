// substitution matrix for multiple sequence alignments

#pragma once

#include "mas.h"

#include <string>
#include <istream>
#include <cassert>

// --------------------------------------------------------------------

template<typename T>
class matrix
{
  public:
	typedef T value_type;
	
					matrix(uint32 m, uint32 n)
						: m_m(m)
						, m_n(n)
					{
						m_data = new value_type[m_m * m_n];
					}
					
					matrix(uint32 m, uint32 n, const T& v)
						: m_m(m)
						, m_n(n)
					{
						m_data = new value_type[m_m * m_n];
						std::fill(m_data, m_data + (m_m * m_n), v);
					}
					
	virtual			~matrix()
					{
						delete [] m_data;
					}
	
	value_type		operator()(uint32 i, uint32 j) const
					{
						assert(i < m_m); assert(j < m_n);
						return m_data[i + j * m_m];
					}
					
	value_type&		operator()(uint32 i, uint32 j)
					{
						assert(i < m_m); assert(j < m_n);
						return m_data[i + j * m_m];
					}

	void			print(std::ostream& os) const
					{
						for (uint32 x = 0; x < m_m; ++x)
						{
							for (uint32 y = 0; y < m_n; ++y)
								os << double(m_data[x + y * m_m]) << ';';
							os << std::endl;
						}
					}

  private:
					matrix(const matrix&);
	matrix&			operator=(const matrix&);

	value_type*		m_data;
	uint32			m_m, m_n;
};

template<typename T>
std::ostream& operator<<(std::ostream& lhs, matrix<T>& rhs)
{
	rhs.print(lhs); return lhs;
}

// --------------------------------------------------------------------

template<typename T>
class symmetric_matrix
{
  public:
	typedef T value_type;
	
					symmetric_matrix(uint32 n);
	virtual			~symmetric_matrix();
	
	value_type		operator()(uint32 i, uint32 j) const;
	value_type&		operator()(uint32 i, uint32 j);
	
	// erase two rows, add one at the end (for neighbour joining)
	void			erase_2(uint32 i, uint32 j);

	void			print(std::ostream& os) const;

  private:
	value_type*		m_data;
	uint32			m_n;
};

template<typename T>
symmetric_matrix<T>::symmetric_matrix(uint32 n)
	: m_n(n)
{
	m_data = new value_type[(m_n * (m_n - 1)) / 2];
}

template<typename T>
symmetric_matrix<T>::~symmetric_matrix()
{
	delete[] m_data;
}

template<typename T>
inline
T& symmetric_matrix<T>::operator()(uint32 i, uint32 j)
{
	if (i > j)
		std::swap(i, j);
	
	assert(j < m_n); assert(i != j);
	return m_data[(j * (j - 1)) / 2 + i];
}

template<typename T>
inline
T symmetric_matrix<T>::operator()(uint32 i, uint32 j) const
{
	if (i > j)
		std::swap(i, j);
	
	assert(j < m_n); assert(i != j);
	return m_data[(j * (j - 1)) / 2 + i];
}

template<typename T>
void symmetric_matrix<T>::erase_2(uint32 di, uint32 dj)
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
void symmetric_matrix<T>::print(std::ostream& os) const
{
	for (uint32 y = 1; y < m_n; ++y)
	{
		os << std::setw(5) << y;
		
		for (uint32 x = 0; x < y; ++x)
			os << (boost::format("  %1.2f") % operator()(x, y));
		
		os << std::endl;
	}
}

template<typename T>
std::ostream& operator<<(std::ostream& lhs, const symmetric_matrix<T>& rhs)
{
	rhs.print(lhs);
	return lhs;
}

// --------------------------------------------------------------------

class substitution_matrix
{
  public:
						substitution_matrix(const std::string& name);
						substitution_matrix(
							const substitution_matrix& m, bool positive);

	virtual				~substitution_matrix() {}

	int8				operator()(aa a, aa b) const
						{
							return m_matrix(a, b);
						}

	float				mismatch_average() const		{ return m_mismatch_average; }
	float				scale_factor() const			{ return m_scale_factor; }

  private:
						substitution_matrix(
							const substitution_matrix&);
	substitution_matrix&
						operator=(const substitution_matrix&);

	void				read(std::istream& is);

	matrix<int8>		m_matrix;
	float				m_mismatch_average;
	float				m_scale_factor;
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
							
							uint32 ix = 0;
							while (distance < m_cutoff[ix] and ix < 3)
								++ix;
							
							if (positive)
								result = m_pos_smat[ix];
							else
								result = m_smat[ix];
							
							return *result;
						}

  private:
						substitution_matrix_family(
							const substitution_matrix_family&);
	substitution_matrix_family&
						operator=(const substitution_matrix_family&);

	float				m_cutoff[4];
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
