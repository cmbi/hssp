#ifndef XSSP_FETCHDBREFS_H
#define XSSP_FETCHDBREFS_H

#pragma once

#include <map>
#include <string>
#include <vector>

void FetchPDBReferences(const std::string& inBaseURL,
  const std::string& inDb, const std::string& inID,
  std::vector<std::string>& outReferences);

void FetchPDBReferences(const std::string& inBaseURL,
  const std::string& inDb,
  std::map<std::string,std::vector<std::string>>& ioReferences);

#endif
