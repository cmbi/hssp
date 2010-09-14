// mas is a reimplementation of Clustal W
//

#pragma once

#if defined(_MSC_VER)
#include <ciso646>
#define snprintf _snprintf
#pragma warning (disable : 4996)
#pragma warning (disable : 4355)
#endif

#include <string>
#include <cassert>
#include <ostream>
#include <iomanip>

#include <boost/format.hpp>

typedef char			int8;
typedef unsigned char	uint8;
typedef short			int16;
typedef unsigned short	uint16;
typedef long			int32;
typedef unsigned long	uint32;

extern int VERBOSE;

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

struct entry
{
					entry(uint32 nr, const std::string& id, const sequence& seq, float weight = 1.0f)
						: m_nr(nr)
						, m_id(id)
						, m_seq(seq)
						, m_weight(weight) {}

	uint32			nr() const						{ return m_nr; }
	float			weight() const					{ return m_weight; }
	uint32			length() const					{ return m_seq.length(); }

	void			insert_gap(uint32 pos);
	void			append_gap();

	void			dump_positions()				{ m_positions.clear(); }

	uint32			m_nr;
	std::string		m_id;
	sequence		m_seq;
	float			m_weight;
	std::vector<uint16>
					m_positions;
};

struct alignment
{
	uint32				length() const				{ return m_entries.front()->m_seq.length(); }

	std::vector<entry*>	m_entries;
};

