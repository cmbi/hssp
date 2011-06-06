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
	MProtein&					inProtein,
	std::ostream&				outHSSP);

}
