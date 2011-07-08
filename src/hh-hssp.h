// maxhom version of hssp generating code
//
//	Copyright, M.L. Hekkelman, UMC St. Radboud, Nijmegen
//

#pragma once

#include <iostream>

class MProtein;

namespace hh
{

void CreateHSSP(
	CDatabankPtr				inDatabank,
	const std::string&			inClustalO,
	MProtein&					inProtein,
	std::ostream&				outHSSP);

void CreateHSSP(
	CDatabankPtr				inDatabank,
	const std::string&			inClustalO,
	const std::string&			inProtein,
	std::ostream&				outHSSP);

}
