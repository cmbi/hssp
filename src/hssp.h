#pragma once

#include <vector>
#include <ostream>

class MProtein;

void CreateHSSPForAlignments(MProtein& inProtein,
	uint32 inMaxHits, uint32 inMinLength, float inCutOff,
	std::vector<std::string>& inChainAlignments, std::ostream& outHSSP);
