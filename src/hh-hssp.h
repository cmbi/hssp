// maxhom version of hssp generating code
//
//	Copyright, M.L. Hekkelman, UMC St. Radboud, Nijmegen
//

#pragma once

#include <iostream>
#include <vector>

class MProtein;

namespace hh
{
	
struct seq
{
	std::string		m_id;
	std::string		m_desc;
	std::string		m_seq;
	
					seq() {}

					seq(const std::string& id, const std::string& seq)
						: m_id(id)
						, m_seq(seq)
					{
					}
	
					seq(const std::string& id, const std::string& seq, const std::string& desc)
						: m_id(id)
						, m_desc(desc)
						, m_seq(seq)
					{
					}
};

typedef std::vector<seq> mseq;

void CreateHSSP(
	CDatabankPtr				inDatabank,
	const std::string&			inClustalO,
	MProtein&					inProtein,
	std::ostream&				outHSSP);

void CreateHSSP(
	CDatabankPtr				inDatabank,
	const std::string&			inClustalO,
	const std::string&			inProtein,
	std::ostream&				outHSSP);

void CreateHSSPForAlignment(
	const mseq&					inAlignment,
	MProtein&					inProtein,
	char						inChain,
	std::ostream&				outHSSP);

}
