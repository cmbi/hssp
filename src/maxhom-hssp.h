// maxhom version of hssp generating code
//
//	Copyright, M.L. Hekkelman, UMC St. Radboud, Nijmegen
//

#pragma once

#include <iostream>

class MProtein;

namespace maxhom
{

void CreateHSSP(
	CDatabankPtr				inDatabank,
	const std::string&			inMaxHom,
	MProtein&					inProtein,
	std::ostream&				outHSSP);

void GetHSSPForHitsAndDSSP(
	CDatabankPtr				inDatabank,
	const std::string&			inMaxHom,
	const std::string&			inPDBID,
	const std::vector<uint32>&	inHits,
	const std::string&			inDSSP,
	int							inMaxAlign,
	std::ostream&				outHSSP);

}
