// Copyright Maarten L. Hekkelman, Radboud University 2008-2011.
//   Distributed under the Boost Software License, Version 1.0.
//       (See accompanying file LICENSE_1_0.txt or copy at    
//             http://www.boost.org/LICENSE_1_0.txt)      
// 
// maxhom version of hssp generating code
//
//	Copyright, M.L. Hekkelman, UMC St. Radboud, Nijmegen
//

#ifndef XSSP_MAXHOM_H
#define XSSP_MAXHOM_H

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

#endif
