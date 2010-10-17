// align-3d

#pragma once

class substitution_matrix_family;

void align_structures(
	const std::string& structureA, const std::string& structureB,
	uint32 iterations,
	substitution_matrix_family& mat, float gop, float gep, float magic);

void test_ss(const std::string& inID);
