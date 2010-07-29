// align is a reimplementation of Clustal W
//

#pragma once

#include <string>
#include <cassert>
#include <ostream>
#include <iomanip>

#include <boost/format.hpp>

#if defined(_MSC_VER)
#include <ciso646>
#define snprintf _snprintf
#endif

typedef char			int8;
typedef unsigned char	uint8;
typedef short			int16;
typedef unsigned short	uint16;
typedef long			int32;
typedef unsigned long	uint32;

#if defined(DEBUG) || ! defined(NDEBUG)
extern int DEBUG;
#endif

// --------------------------------------------------------------------

typedef uint8					aa;
typedef std::basic_string<aa>	sequence;

std::string decode(const sequence& s);
sequence encode(const std::string& s);

const uint8 kAA[] = {
	'A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G',
	'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S',
	'T', 'W', 'Y', 'V', 'B', 'Z', 'X', '*',
	'-'
};

const uint32
	kAACount = sizeof(kAA),
	kFilteredCode = 22,
	kUnknownCode = 23,
	kSignalGapCode = 24,
	kSentinalScore = kSignalGapCode;

extern aa kAAReverse[256];

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
						for (uint32 y = 0; y < m_m; ++y)
						{
							os << std::setw(3) << y;
							for (uint32 x = 0; x < m_n; ++x)
								os << ' ' << std::setw(5) << int32(m_data[x + y * m_m]);
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

class my_bad : public std::exception
{
  public:
					my_bad(const std::string& msg);
					my_bad(const boost::format& msg);

	virtual const char*
					what() const throw()	{ return m_msg; }

  private:
	char			m_msg[1024];
};

