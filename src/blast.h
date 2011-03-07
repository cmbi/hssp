// Blast code for hssp/dssp
//
//	Copyright, M.L. Hekkelman, UMC St. Radboud, Nijmegen
//

#pragma once

#include "CDatabank.h"

extern int BLAST_THREADS;

void BlastSequence(
	CDatabankPtr			inDatabank,
	const std::string&		inSequence,
	std::vector<uint32>&	outHits);
