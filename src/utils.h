// utility routines for mas

#pragma once

#ifndef NDEBUG
#include <iostream>
#endif

#include <time.h>
#include <boost/thread.hpp>

// --------------------------------------------------------------------

class mas_exception : public std::exception
{
  public:
					mas_exception(const std::string& msg);
					mas_exception(const boost::format& msg);

	virtual const char*
					what() const throw()	{ return m_msg; }

  private:
	char			m_msg[1024];
};

// --------------------------------------------------------------------

class progress
{
  public:
					progress(const std::string& message, uint32 max);
					~progress();

	void			step(uint32 advance = 1);

  private:
	
	void			run();

	std::string		m_msg;
	uint32			m_step, m_max;
	boost::mutex	m_mutex;
	boost::thread	m_thread;
	std::time_t		m_start;
};

// --------------------------------------------------------------------

#ifndef NDEBUG
struct stats
{
	stats() : m_max(0), m_count(0), m_cumm(0) {}
	~stats();
	
	void operator()(uint32 i)
	{
		if (m_max < i)
			m_max = i;
		++m_count;
		m_cumm += i;
	}
	
	uint32 m_max, m_count, m_cumm;
};
#endif
