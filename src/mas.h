// Copyright Maarten L. Hekkelman, Radboud University 2008-2011.
//   Distributed under the Boost Software License, Version 1.0.
//       (See accompanying file LICENSE_1_0.txt or copy at    
//             http://www.boost.org/LICENSE_1_0.txt)      

#pragma once

#include <string>

#if defined(_MSC_VER)

#ifndef P_WIN
#define P_WIN 1
#endif

// These are Microsoft Visual C++ special settings
// the iso646 file contains the C++ keywords that are
// otherwise not recognized.
#include <ciso646>
#define snprintf _snprintf

// Disable some warnings
#pragma warning (disable : 4996)
#pragma warning (disable : 4355)
#endif

#include <boost/version.hpp>
#include <boost/cstdint.hpp>

typedef int8_t		int8;
typedef uint8_t		uint8;
typedef int16_t		int16;
typedef uint16_t	uint16;
typedef int32_t		int32;
typedef uint32_t	uint32;
typedef int64_t		int64;
typedef uint64_t	uint64;

#ifndef nullptr
#define nullptr NULL
#endif

// we even have globals:
extern int VERBOSE;

// Code for amino acid sequences

typedef std::basic_string<uint8> sequence;

// 22 real letters and 1 dummy (X is the dummy, B and Z are pseudo letters)
extern const char kResidues[]; // = "ACDEFGHIKLMNPQRSTVWYBZX";
extern const uint8 kResidueNrTable[];

inline uint8 ResidueNr(char inAA)
{
	int result = 23;

	inAA |= 040;
	if (inAA >= 'a' and inAA <= 'z')
		result = kResidueNrTable[inAA - 'a'];

	return result;
}

inline bool is_gap(char aa)
{
	return aa == ' ' or aa == '.' or aa == '-';
}

sequence encode(const std::string& s);
std::string decode(const sequence& s);
