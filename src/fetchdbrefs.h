#pragma once

#include <vector>
#include <string>

void FetchPDBReferences(const std::string& inBaseURL,
	const std::string& inDb, const std::string& inID,
	std::vector<std::string>& outReferences);
