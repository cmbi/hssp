#pragma once

#include <vector>

#include <boost/filesystem/path.hpp>

class MProtein;

namespace HSSP
{
	
extern const float kThreshold, kFragmentCutOff;

void CreateHSSP(const MProtein& inProtein,
	const std::vector<boost::filesystem::path>& inDatabanks,
	uint32 inMaxhits, uint32 inMinSeqLength, float inGapOpen, float inGapExtend,
	float inThreshold, float inFragmentCutOff, uint32 inThreads,
	std::ostream& inOutStream);

void CreateHSSP(const std::string& inProtein,
	const std::vector<boost::filesystem::path>& inDatabanks,
	uint32 inMaxhits, uint32 inMinSeqLength, float inGapOpen, float inGapExtend,
	float inThreshold, float inFragmentCutOff, uint32 inThreads,
	std::ostream& inOutStream);

}
