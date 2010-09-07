// i/o code for mas

#pragma once

#include <boost/filesystem/path.hpp>
#include <vector>
#include <string>

struct entry;

void readFasta(boost::filesystem::path path, std::vector<entry>& seq);
void readAlignmentFromHsspFile(boost::filesystem::path path,
	char& chainID, std::vector<entry>& seq);

void report(const std::vector<entry*>& alignment,
	std::ostream& os, const std::string& format);
