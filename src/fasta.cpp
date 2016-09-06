#include "fasta.h"
#include "structure.h"
#include "utils.h"

#include <boost/algorithm/string.hpp>

namespace ba = boost::algorithm;


// TODO: use unique_ptr
std::vector<MProtein*> read_proteins_from_fasta(std::istream& in)
{
  std::vector<MProtein*> proteins;
  for (std::string line; std::getline(in, line);)
  {
    // Check for identifier line
    if (!ba::starts_with(line, ">"))
      throw mas_exception("FASTA id line does not start with '>'");

    // Parse the identifier
    std::string id = line;
    id.erase(0, 1);
    std::string::size_type s = id.find(' ');
    if (s != std::string::npos)
      id.erase(s);

    if (id.empty())
      throw mas_exception("Missing FASTA id");

    // Parse the data for current protein
    std::streambuf* b = in.rdbuf();
    std::string seq;
    while (b->sgetc() != std::streambuf::traits_type::eof() &&
           b->sgetc() != '>')
    {
      std::string line;
      std::getline(in, line);
      seq += line;
    }

    // Create the Protein instance
    MChain* chain = new MChain("A");
    std::vector<MResidue*>& residues = chain->GetResidues();
    MResidue* last = nullptr;
    int32 nr = 1;
    for (auto& r: seq)
    {
      residues.push_back(new MResidue(nr, r, last));
      ++nr;
      last = residues.back();
    }

    MProtein* result = new MProtein(id, chain);
    proteins.push_back(result);
  }

  return proteins;
}
