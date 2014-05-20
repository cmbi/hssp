// Copyright Maarten L. Hekkelman, Radboud University 2008-2011.
//   Distributed under the Boost Software License, Version 1.0.
//       (See accompanying file LICENSE_1_0.txt or copy at    
//             http://www.boost.org/LICENSE_1_0.txt)      
// 
// mas guide tree support routines and class

#ifndef XSSP_GUIDE_H
#define XSSP_GUIDE_H

#pragma once

#include <iostream>
#include <map>
#include <boost/tr1/tuple.hpp>

#include "align-2d.h"

// --------------------------------------------------------------------

class GuideTreeParser
{
  public:
						GuideTreeParser(std::istream& data, std::map<std::string,leaf_node*>& m)
							: m_data(data)
							, m_map(m)
						{
							getNextToken();
						}
	
	base_node*			parse();

  private:

	enum GuideTreeToken {
		gtt_Undef,
		gtt_Open,
		gtt_Close,
		gtt_End,
		gtt_ID,
		gtt_Colon,
		gtt_Comma,
		gtt_Weight
	};
	
	void				getNextToken();
	void				match(GuideTreeToken token);
	
	std::tr1::tuple<base_node*,float>
						parseNode();
	base_node*			parseGroup();
	
	std::istream&		m_data;
	std::map<std::string,leaf_node*>&
						m_map;
	float				m_value;
	std::string			m_token;
	GuideTreeToken		m_lookahead;
};

void useGuideTree(const std::string& guide, std::vector<base_node*>& tree);

#endif
