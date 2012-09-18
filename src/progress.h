// Copyright Maarten L. Hekkelman, Radboud University 2008-2011.
//   Distributed under the Boost Software License, Version 1.0.
//       (See accompanying file LICENSE_1_0.txt or copy at    
//             http://www.boost.org/LICENSE_1_0.txt)      

#pragma once

#if defined(__linux__)

#if defined(__INTEL_COMPILER_BUILD_DATE) || (defined(__GNUC__) && (__GNUC__ > 4 || ( __GNUC__ == 4 && __GNUC_MINOR__ >= 6)))
#include <atomic>
typedef std::atomic<int64>	MCounter;
#else
#include <boost/detail/atomic_count.hpp>
typedef boost::detail::atomic_count        MCounter;
#endif

#else

struct MCounter
{
	MCounter(int64 inValue) : m_value(inValue) {}

			operator int64() const					{ return m_value; }

	int64	operator++(int);
	int64	operator+=(int64 inValue);
	int64	operator=(int64 inValue);
//	bool	operator==(const MCounter& rhs) const 	{ return m_value == rhs.m_value; }

	int64	m_value;
};


//typedef int64 MCounter;
//
//int64 add(MCounter& ioCounter, int64 inIncrement);
//int64 set(MCounter& ioCounter, int64 inValue);

#endif

// --------------------------------------------------------------------

class MProgress
{
  public:
				MProgress(int64 inMax, const std::string& inAction);
	virtual		~MProgress();
	
	void		Consumed(int64 inConsumed);	// consumed is relative
	void		Progress(int64 inProgress);	// progress is absolute

	void		Message(const std::string& inMessage);

  private:
				MProgress(const MProgress&);
	MProgress&	operator=(const MProgress&);

	struct MProgressImpl*	mImpl;
};
