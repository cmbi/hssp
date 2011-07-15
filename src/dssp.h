// A DSSP reimplementation
//
//	Copyright, M.L. Hekkelman, UMC St. Radboud, Nijmegen
//

#pragma once

#include <iostream>

class MProtein;
class MResidue;

std::string ResidueToDSSPLine(const MResidue& residue);
void WriteDSSP(MProtein& protein, std::ostream& os);
