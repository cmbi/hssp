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
#include <boost/cstdint.hpp>

typedef int8_t		int8;
typedef uint8_t		uint8;
typedef int16_t		int16;
typedef uint16_t	uint16;
typedef int32_t		int32;
typedef uint32_t	uint32;
typedef int64_t		int64;
typedef uint64_t	uint64;

#ifndef nil
#define nil NULL
#endif

extern int VERBOSE, MULTI_THREADED;

