// maxhom version of hssp generating code
//
//	Copyright, M.L. Hekkelman, UMC St. Radboud, Nijmegen
//

#pragma once

#include <iostream>
#include <vector>

class MProtein;

namespace hmmer
{
	
void CreateHSSP(
	CDatabankPtr				inDatabank,
	MProtein&					inProtein,
	const std::string&			inJackHmmer,
	uint32						inIterations,
	std::ostream&				outHSSP);

void CreateHSSP(
	CDatabankPtr				inDatabank,
	const std::string&			inProtein,
	const std::string&			inJackHmmer,
	uint32						inIterations,
	std::ostream&				outHSSP);

//void CreateHSSP(
//	CDatabankPtr				inDatabank,
//	const std::string&			inClustalO,
//	MProtein&					inProtein,
//	std::ostream&				outHSSP);
//
//void CreateHSSP(
//	CDatabankPtr				inDatabank,
//	const std::string&			inClustalO,
//	const std::string&			inProtein,
//	std::ostream&				outHSSP);


}
