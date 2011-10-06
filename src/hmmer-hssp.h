// version of hssp using jackhmmer
//
//	Copyright, M.L. Hekkelman, UMC St. Radboud, Nijmegen
//

#pragma once

#include <iostream>
#include <vector>

#include <boost/filesystem/path.hpp>

#include "CDatabank.h"

class MProtein;

namespace hmmer
{

void CreateHSSP(
	CDatabankPtr					inDatabank,
	MProtein&						inProtein,
	const boost::filesystem::path&	inFastaDir,
	const boost::filesystem::path&	inJackHmmer,
	uint32							inIterations,
	uint32							inMaxHits,
	uint32							inMinSeqLength,
	std::ostream&					outHSSP);

void CreateHSSP(
	CDatabankPtr					inDatabank,
	const std::string&				inProtein,
	const boost::filesystem::path&	inFastaDir,
	const boost::filesystem::path&	inJackHmmer,
	uint32							inIterations,
	uint32							inMaxHits,
	std::ostream&					outHSSP);

void CreateHSSP(
	CDatabankPtr					inDatabank,
	const MProtein&					inProtein,
	const boost::filesystem::path&	inDataDir,		// for the stockholm files
	const boost::filesystem::path&	inFastaDir,
	const boost::filesystem::path&	inJackHmmer,
	uint32							inIterations,
	uint32							inMaxHits,
	std::vector<std::string>		inStockholmIds,
	std::ostream&					outHSSP);

// Create a FastA formatted aligment from a Stockholm file.
// The stockholm file should have been created with the
// inQuery sequence. The alignment is cut to match the query
// (leading and trailing alignment info is discarded).
// The resulting alignment is ordered by HSSP score and
// written out as a FastA file.
// (query sequence may be left empty)
void ConvertHmmerAlignment(
	const std::string&				inQuerySequence,
	const boost::filesystem::path&	inStockholmFile,
	const boost::filesystem::path&	inFastaFile);
	
}
