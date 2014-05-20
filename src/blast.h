// Copyright Maarten L. Hekkelman, Radboud University 2008-2011.
//   Distributed under the Boost Software License, Version 1.0.
//       (See accompanying file LICENSE_1_0.txt or copy at    
//             http://www.boost.org/LICENSE_1_0.txt)      
// 
//	Simplified Blast algorithm implementation. Works on
//	FastA formatted files containing proteins.

#ifndef XSSP_BLAST_H
#define XSSP_BLAST_H

#pragma once

#include <vector>

#include <boost/filesystem/path.hpp>

//struct BlastHit
//{
//	std::string	id;
//	std::string	acc;
//	std::string	def;
//	float		score;
//	sequence	seq;
//};
//
//// blast with default parameters
//void Blast(const std::string& seq,
//	const std::vector<boost::filesystem::path>& db,
//	uint32 inReportLimit, std::vector<BlastHit>& outHits);

void SearchAndWriteResultsAsFastA(std::ostream& inOutFile,
	const std::vector<boost::filesystem::path>& inDatabanks,
	const std::string& inQuery, const std::string& inProgram,
	const std::string& inMatrix, uint32 inWordSize, double inExpect,
	bool inFilter, bool inGapped, int32 inGapOpen, int32 inGapExtend,
	uint32 inReportLimit, uint32 inThreads);

#endif
