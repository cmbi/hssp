// align-3d

#pragma once

class substitution_matrix_family;

void align_structures(
	std::istream& structureA, std::istream& structureB,
	char chainA, char chainB,
	uint32 iterations,
	substitution_matrix_family& mat, float gop, float gep, float magic,
	std::vector<entry*>& outAlignment);
