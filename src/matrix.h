// Copyright Maarten L. Hekkelman, Radboud University 2008-2011.
//   Distributed under the Boost Software License, Version 1.0.
//       (See accompanying file LICENSE_1_0.txt or copy at    
//             http://www.boost.org/LICENSE_1_0.txt)      
// 
// substitution matrix for multiple sequence alignments

#pragma once

#include "mas.h"
#include "align-2d.h"

#include <string>
#include <istream>
#include <cassert>
#include <stdexcept>

// --------------------------------------------------------------------
// uBlas compatible matrix types
// matrix is m x n, addressing i,j is 0 <= i < m and 0 <= j < n
// element m i,j is mapped to [i * n + j] and thus storage is row major

template<typename T>
class matrix_base
{
  public:

	typedef T value_type;

	virtual				~matrix_base() {}

	virtual uint32		dim_m() const = 0;
	virtual uint32		dim_n() const = 0;

	virtual value_type&	operator()(uint32 i, uint32 j) { throw std::runtime_error("unimplemented method"); }
	virtual value_type	operator()(uint32 i, uint32 j) const = 0;
	
	matrix_base&		operator*=(const value_type& rhs);

	matrix_base&		operator-=(const value_type& rhs);
};

template<typename T>
matrix_base<T>& matrix_base<T>::operator*=(const T& rhs)
{
	for (uint32 i = 0; i < dim_m(); ++i)
	{
		for (uint32 j = 0; j < dim_n(); ++j)
		{
			operator()(i, j) *= rhs;
		}
	}
	
	return *this;
}

template<typename T>
matrix_base<T>& matrix_base<T>::operator-=(const T& rhs)
{
	for (uint32 i = 0; i < dim_m(); ++i)
	{
		for (uint32 j = 0; j < dim_n(); ++j)
		{
			operator()(i, j) -= rhs;
		}
	}
	
	return *this;
}

template<typename T>
std::ostream& operator<<(std::ostream& lhs, const matrix_base<T>& rhs)
{
	lhs << '[' << rhs.dim_m() << ',' << rhs.dim_n() << ']' << '(';
	for (uint32 i = 0; i < rhs.dim_m(); ++i)
	{
		lhs << '(';
		for (uint32 j = 0; j < rhs.dim_n(); ++j)
		{
			if (j > 0)
				lhs << ',';
			lhs << rhs(i,j);
		}
		lhs << ')';
	}
	lhs << ')';
	
	return lhs;
}

template<typename T>
class matrix : public matrix_base<T>
{
  public:
	typedef T value_type;

						template<typename T2>
						matrix(const matrix_base<T2>& m)
							: m_m(m.dim_m())
							, m_n(m.dim_n())
						{
							m_data = new value_type[m_m * m_n];
							for (uint32 i = 0; i < m_m; ++i)
							{
								for (uint32 j = 0; j < m_n; ++j)
									operator()(i, j) = m(i, j);
							}
						}

						matrix(const matrix& m)
							: m_m(m.m_m)
							, m_n(m.m_n)
						{
							m_data = new value_type[m_m * m_n];
							std::copy(m.m_data, m.m_data + (m_m * m_n), m_data);
						}

	matrix&				operator=(const matrix& m)
						{
							value_type t = new value_type[m.m_m * m.m_n];
							std::copy(m.m_data, m.m_data + (m_m * m_n), t);
							
							delete[] m_data;
							m_data = t;
							m_m = m.m_m;
							m_n = m.m_n;
							
							return *this;
						}
	
						matrix(uint32 m, uint32 n, T v = T())
							: m_m(m)
							, m_n(n)
						{
							m_data = new value_type[m_m * m_n];
							std::fill(m_data, m_data + (m_m * m_n), v);
						}
						
	virtual				~matrix()
						{
							delete [] m_data;
						}
	
	virtual uint32		dim_m() const 					{ return m_m; }
	virtual uint32		dim_n() const					{ return m_n; }

	virtual value_type	operator()(uint32 i, uint32 j) const
						{
							assert(i < m_m); assert(j < m_n);
							return m_data[i * m_n + j];
						}
					
	virtual value_type&	operator()(uint32 i, uint32 j)
						{
							assert(i < m_m); assert(j < m_n);
							return m_data[i * m_n + j];
						}

  private:
	value_type*			m_data;
	uint32				m_m, m_n;
};

// --------------------------------------------------------------------

template<typename T>
class symmetric_matrix : public matrix_base<T>
{
  public:
	typedef typename matrix_base<T>::value_type value_type;

						symmetric_matrix(uint32 n)
							: m_owner(true)
							, m_n(n)
						{
							uint32 N = (m_n * (m_n + 1)) / 2;
							m_data = new value_type[N];
							std::fill(m_data, m_data + N, T(0));
						}

						symmetric_matrix(const T* data, uint32 n)
							: m_owner(false)
							, m_data(const_cast<T*>(data))
							, m_n(n)
						{
						}


	virtual				~symmetric_matrix()
						{
							if (m_owner)
								delete[] m_data;
						}
	
	virtual uint32		dim_m() const					{ return m_n; }
	virtual uint32		dim_n() const					{ return m_n; }

	T					operator()(uint32 i, uint32 j) const;
	virtual T&			operator()(uint32 i, uint32 j);
	
	// erase two rows, add one at the end (for neighbour joining)
	void				erase_2(uint32 i, uint32 j);

  private:
	bool				m_owner;
	value_type*			m_data;
	uint32				m_n;
};

template<typename T>
inline
T symmetric_matrix<T>::operator()(uint32 i, uint32 j) const
{
	return i < j
		? m_data[(j * (j + 1)) / 2 + i]
		: m_data[(i * (i + 1)) / 2 + j];
//	if (i > j)
//		std::swap(i, j);
//	assert(j < m_n);
//	return m_data[(j * (j + 1)) / 2 + i];
}

template<typename T>
inline
T& symmetric_matrix<T>::operator()(uint32 i, uint32 j)
{
	if (i > j)
		std::swap(i, j);
	assert(j < m_n);
	return m_data[(j * (j + 1)) / 2 + i];
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
class identity_matrix : public matrix_base<T>
{
  public:
	typedef typename matrix_base<T>::value_type value_type;

						identity_matrix(uint32 n)
							: m_n(n)
						{
						}

	virtual uint32		dim_m() const					{ return m_n; }
	virtual uint32		dim_n() const					{ return m_n; }

	virtual value_type	operator()(uint32 i, uint32 j) const
						{
							value_type result = 0;
							if (i == j)
								result = 1;
							return result;
						}

  private:
	uint32				m_n;
};

// --------------------------------------------------------------------
// matrix functions

template<typename T>
matrix<T> operator*(const matrix_base<T>& lhs, const matrix_base<T>& rhs)
{
	matrix<T> result(min(lhs.dim_m(), rhs.dim_m()), min(lhs.dim_n(), rhs.dim_n()));
	
	for (uint32 i = 0; i < result.dim_m(); ++i)
	{
		for (uint32 j = 0; j < result.dim_n(); ++j)
		{
			for (uint32 li = 0, rj = 0; li < lhs.dim_m() and rj < rhs.dim_n(); ++li, ++rj)
				result(i, j) += lhs(li, j) * rhs(i, rj);
		}
	}
	
	return result;
}

template<typename T>
matrix<T> operator*(const matrix_base<T>& lhs, T rhs)
{
	matrix<T> result(lhs);
	result *= rhs;

	return result;
}

template<typename T>
matrix<T> operator-(const matrix_base<T>& lhs, const matrix_base<T>& rhs)
{
	matrix<T> result(min(lhs.dim_m(), rhs.dim_m()), min(lhs.dim_n(), rhs.dim_n()));
	
	for (uint32 i = 0; i < result.dim_m(); ++i)
	{
		for (uint32 j = 0; j < result.dim_n(); ++j)
		{
			result(i, j) = lhs(i, j) - rhs(i, j);
		}
	}
	
	return result;
}

template<typename T>
matrix<T> operator-(const matrix_base<T>& lhs, T rhs)
{
	matrix<T> result(lhs.dim_m(), lhs.dim_n());
	result -= rhs;
	return result;
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

	int8				operator()(char a, char b) const
						{
							return m_matrix(encode(a), encode(b));
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
