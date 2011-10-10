// Copyright Maarten L. Hekkelman, Radboud University 2008-2011.
//   Distributed under the Boost Software License, Version 1.0.
//       (See accompanying file LICENSE_1_0.txt or copy at    
//             http://www.boost.org/LICENSE_1_0.txt)      
//
// buffer is a thread safe queue

#pragma once

#include <deque>
#include <boost/thread/mutex.hpp>
#include <boost/thread/condition.hpp>

template<class T, uint32 N = 100>
class buffer
{
  public:

						buffer() {}

	void				put(T inValue);
	T					get();

  private:
						buffer(const buffer&);
	buffer&				operator=(const buffer&);

	std::deque<T>		m_queue;
	boost::mutex		m_mutex;
	boost::condition	m_empty, m_full;
};

template<class T, uint32 N>
void buffer<T,N>::put(T inValue)
{
	boost::mutex::scoped_lock lock(m_mutex);

	while (m_queue.size() >= N)
		m_full.wait(lock);
	
	m_queue.push_back(inValue);

	m_empty.notify_one();
}

template<class T, uint32 N>
T buffer<T,N>::get()
{
	boost::mutex::scoped_lock lock(m_mutex);

	while (m_queue.empty())
		m_empty.wait(lock);
	
	T result = m_queue.front();
	m_queue.pop_front();

	m_full.notify_one();
	
	return result;
}
