// utility routines for mas

#include "mas.h"

#include <iomanip>
#include <iostream>
#include <boost/thread.hpp>
#include <time.h>
#include <termios.h>
#include <sys/ioctl.h>

#include "utils.h"

using namespace std;

// --------------------------------------------------------------------

progress::progress(const string& msg, uint32 max)
	: m_msg(msg)
	, m_step(0)
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

			// first, get the current terminal width
		    struct winsize w;
		    ioctl(0, TIOCGWINSZ, &w);

			string msg;
			if (m_msg.length() > w.ws_col)
				msg = m_msg.substr(0, w.ws_col);
			else
			{
				msg = m_msg;
				
				if (msg.length() + 10 < w.ws_col)
				{
					float fraction = float(m_step) / m_max;
					if (fraction < 0)
						fraction = 0;
					else if (fraction > 1)
						fraction = 1;
					
					uint32 thermometer_width = w.ws_col - msg.length() - 10;
					uint32 thermometer_done_width = uint32(fraction * thermometer_width);
					msg += " [";
					msg += string(thermometer_done_width, '#');
					msg += string(thermometer_width - thermometer_done_width, '-');
					msg += ']';
					msg += string(w.ws_col - msg.length(), ' ');
				}
			}

			cout << '\r' << msg << flush;
			first = false;
		}
	}
	catch (boost::thread_interrupted& e)
	{
	    struct winsize w;
	    ioctl(0, TIOCGWINSZ, &w);

		if (not first)
			cout << '\r' << m_msg << " done" << string(w.ws_col - m_msg.length() - 5, ' ') << endl;
	}
}
