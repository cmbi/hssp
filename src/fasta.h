#ifndef XSSP_FASTA_H
#define XSSP_FASTA_H

#pragma once

#include <istream>
#include <vector>

class MProtein;


std::vector<MProtein*> read_proteins_from_fasta(std::istream& in);


#endif
