// utility routines for mas

#include "align.h"

#include <iomanip>

#include "utils.h"
#include <boost/thread.hpp>
#include <time.h>

// --------------------------------------------------------------------

progress::progress(uint32 max)
	: m_step(0)
	, m_max(max)
	, m_thread(boost::bind(&progress::run, this))
{
	time(&m_start);
}

progress::~progress()
{
	boost::mutex::scoped_lock lock(m_mutex);
	m_step = m_max;
	m_thread.interrupt();
	m_thread.join();
}

void progress::step()
{
	boost::mutex::scoped_lock lock(m_mutex);
	++m_step;
}

void progress::run()
{
	bool first = true;

	try
	{
		for (;;)
		{
			boost::this_thread::sleep(boost::posix_time::seconds(1));
			cout << '\r' << m_msg << ' ' << m_step << " of " << m_max << flush;
			first = false;
		}
	}
	catch (boost::thread_interrupted& e)
	{
		if (not first)
			cout << '\r' << "done" << string(16, ' ') << endl;
	}
}
