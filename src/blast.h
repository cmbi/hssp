// Copyright Maarten L. Hekkelman, Radboud University 2008-2011.
//   Distributed under the Boost Software License, Version 1.0.
//       (See accompanying file LICENSE_1_0.txt or copy at    
//             http://www.boost.org/LICENSE_1_0.txt)      
// 
// Blast code for hssp/dssp
//
//	Copyright, M.L. Hekkelman, UMC St. Radboud, Nijmegen
//

#pragma once

#include "CDatabank.h"

class MProtein;

extern int BLAST_THREADS;

void BlastSequence(
	CDatabankPtr			inDatabank,
	const std::string&		inSequence,
	std::vector<uint32>&	outHits);

void BlastProtein(
	CDatabankPtr			inDatabank,
	const MProtein&			inProtein,
	std::vector<uint32>&	outHits);
