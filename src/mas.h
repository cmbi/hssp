// mas is a reimplementation of Clustal W with support for
// predefined blocks of aligned positions in the input sequence.

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

typedef char				int8;
typedef unsigned char		uint8;
typedef short				int16;
typedef unsigned short		uint16;
typedef int					int32;
typedef unsigned int		uint32;
typedef long long			int64;
typedef unsigned long long	uint64;

#ifndef nil
#define nil NULL
#endif

extern int VERBOSE, MULTI_THREADED;

