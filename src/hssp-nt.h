#ifndef HSSP_HSSPNT_H
#define HSSP_HSSPNT_H

#pragma once

#include "mas.h"

#include <boost/filesystem/path.hpp>

#include <vector>


class MProtein;

namespace HSSP
{

extern const float kThreshold, kFragmentCutOff;

void CreateHSSP(const MProtein& inProtein,
  const std::vector<boost::filesystem::path>& inDatabanks,
  uint32 inMaxhits, uint32 inMinSeqLength, float inGapOpen, float inGapExtend,
  float inThreshold, float inFragmentCutOff, uint32 inThreads,
  bool inFetchDBRefs, std::ostream& inOutStream);

void CreateHSSP(const std::string& inProtein,
  const std::vector<boost::filesystem::path>& inDatabanks,
  uint32 inMaxhits, uint32 inMinSeqLength, float inGapOpen, float inGapExtend,
  float inThreshold, float inFragmentCutOff, uint32 inThreads,
  bool inFetchDBRefs, std::ostream& inOutStream);

}

#endif
