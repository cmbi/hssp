#pragma once

#include <vector>
#include <string>

void FetchPDBReferences(const std::string& inBaseURL,
	std::vector<std::string>& outReferences);
