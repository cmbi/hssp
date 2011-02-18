// A DSSP reimplementation
//
//	Copyright, M.L. Hekkelman, UMC St. Radboud, Nijmegen
//

#pragma once

#include <iostream>

class MProtein;

void WriteDSSP(MProtein& protein, std::ostream& os);
