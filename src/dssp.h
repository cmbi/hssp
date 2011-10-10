// Copyright Maarten L. Hekkelman, Radboud University 2008-2011.
//   Distributed under the Boost Software License, Version 1.0.
//       (See accompanying file LICENSE_1_0.txt or copy at    
//             http://www.boost.org/LICENSE_1_0.txt)      

#pragma once

#include <iostream>

class MProtein;
class MResidue;

// Write the DSSP line for a single residue
std::string ResidueToDSSPLine(const MResidue& residue);

// Write a complete DSSP file for a protein
void WriteDSSP(MProtein& protein, std::ostream& os);
