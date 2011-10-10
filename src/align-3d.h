// Copyright Maarten L. Hekkelman, Radboud University 2008-2011.
//   Distributed under the Boost Software License, Version 1.0.
//       (See accompanying file LICENSE_1_0.txt or copy at    
//             http://www.boost.org/LICENSE_1_0.txt)      
// 
// align-3d

#pragma once

class substitution_matrix_family;

void align_structures(
	std::istream& structureA, std::istream& structureB,
	char chainA, char chainB,
	uint32 iterations,
	substitution_matrix_family& mat, float gop, float gep, float magic,
	std::vector<entry*>& outAlignment);
