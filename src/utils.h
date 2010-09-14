// utility routines for mas

#pragma once

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

