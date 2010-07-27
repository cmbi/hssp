// align is a reimplementation of Clustal W
//

#pragma once

#include <string>
#include <cassert>

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
						std::assert(i < m_m); assert(j < m_n);
						return m_data[i + j * m_m];
					}
					
	value_type&		operator()(uint32 i, uint32 j)
					{
						std::assert(i < m_m); assert(j < m_n);
						return m_data[i + j * m_m];
					}

  private:
					matrix(const matrix&);
	matrix&			operator=(const matrix&);

	value_type*		m_data;
	uint32			m_m, m_n;
};

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

