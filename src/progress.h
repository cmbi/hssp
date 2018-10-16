// Copyright Maarten L. Hekkelman, Radboud University 2008-2011.
// Copyright Coos Baakman, Jon Black, Wouter G. Touw & Gert Vriend, Radboud university medical center 2015.
//   Distributed under the Boost Software License, Version 1.0.
//       (See accompanying file LICENSE_1_0.txt or copy at
//             http://www.boost.org/LICENSE_1_0.txt)

#pragma once

#include <atomic>
typedef std::atomic<int64>  MCounter;

// --------------------------------------------------------------------

class MProgress
{
  public:
        MProgress(int64 inMax, const std::string& inAction);
  virtual    ~MProgress();

  void    Consumed(int64 inConsumed);  // consumed is relative
  void    Progress(int64 inProgress);  // progress is absolute

  void    Message(const std::string& inMessage);

  private:
        MProgress(const MProgress&);
  MProgress&  operator=(const MProgress&);

  struct MProgressImpl*  mImpl;
};
