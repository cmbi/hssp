// Copyright Maarten L. Hekkelman, Radboud University 2008-2011.
// Copyright Coos Baakman, Jon Black, Wouter G. Touw & Gert Vriend, Radboud university medical center 2015.
//   Distributed under the Boost Software License, Version 1.0.
//       (See accompanying file LICENSE_1_0.txt or copy at
//             http://www.boost.org/LICENSE_1_0.txt)

#ifndef XSSP_UTILS_H
#define XSSP_UTILS_H

#pragma once

#include "mas.h"

#include <filesystem>
#include <string>

#ifndef NDEBUG
#include <iostream>
#endif
#include <ctime>

// --------------------------------------------------------------------

class arg_vector
{
  public:

        arg_vector(const std::string& program)
        {
          m_args.push_back(program);
        }

  void    push(const std::string& option)
        {
          m_args.push_back(option);
        }

  template<class T>
  void    push(const std::string& option, const T& value);

        operator char* const*();

  private:
  friend std::ostream& operator<<(std::ostream& os, const arg_vector& argv);

  std::vector<std::string>  m_args;
  std::vector<const char*>  m_argv;
};

template<class T>
inline
void arg_vector::push(const std::string& option, const T& value)
{
  m_args.push_back(option);
  m_args.push_back(std::to_string(value));
}

template<>
inline
void arg_vector::push(const std::string& option, const std::string& value)
{
  m_args.push_back(option);
  m_args.push_back(value);
}

std::ostream& operator<<(std::ostream& os, const arg_vector& argv);

// --------------------------------------------------------------------

class mas_exception : public std::exception
{
  public:
          mas_exception(const std::string& msg);

  virtual const char*
          what() const throw()  { return m_msg; }

  private:
  char      m_msg[1024];
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

// --------------------------------------------------------------------

void WriteToFD(int inFD, const std::string& inText);
std::filesystem::path get_home();

std::string Format(const std::string &fmt, ...);

std::string TrimLeft(const std::string &s);
std::string TrimRight(const std::string &s);
std::string Trim(const std::string &s);

bool StartsWith(const std::string &s, const std::string &prefix);


#endif
