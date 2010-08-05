// utility routines for mas

#pragma once

// --------------------------------------------------------------------

class progress
{
  public:
					progress(const std::string& message, uint32 max);
					~progress();

	void			step();

  private:
	
	void			run();

	std::string		m_msg;
	uint32			m_step, m_max;
	boost::mutex	m_mutex;
	boost::thread	m_thread;
	std::time_t		m_start;
};

