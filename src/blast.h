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
