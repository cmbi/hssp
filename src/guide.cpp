// Copyright Maarten L. Hekkelman, Radboud University 2008-2011.
//   Distributed under the Boost Software License, Version 1.0.
//       (See accompanying file LICENSE_1_0.txt or copy at
//             http://www.boost.org/LICENSE_1_0.txt)
//
// guide.cpp - support for DND files / guide trees

#include "mas.h"

#include <iostream>

#include <boost/lexical_cast.hpp>
#include <boost/filesystem/fstream.hpp>
#include <boost/foreach.hpp>
#define foreach BOOST_FOREACH

#include "utils.h"
#include "guide.h"

using namespace std;
using namespace tr1;
namespace fs = boost::filesystem;

void GuideTreeParser::getNextToken()
{
  m_lookahead = gtt_Undef;
  m_token.clear();
  m_value = 0;

  enum State {
    st_Start,
    st_ID,
    st_Number,
    st_Fraction
  } state = st_Start;

  while (m_lookahead == gtt_Undef)
  {
    char ch = 0;
    m_data.get(ch);

    m_token += ch;

    switch (state)
    {
      case st_Start:
        switch (ch)
        {
          case ' ':
          case '\r':
          case '\n':
          case '\t':  m_token.clear();      break;
          case '(':  m_lookahead = gtt_Open;    break;
          case ')':  m_lookahead = gtt_Close;  break;
          case ':':  m_lookahead = gtt_Colon;  break;
          case ',':  m_lookahead = gtt_Comma;  break;
          case ';':  m_lookahead = gtt_End;    break;
          default:
            if (isdigit(ch) or ch == '-')
              state = st_Number;
            else if (isalnum(ch) or ch == '_')
              state = st_ID;
            else
              throw mas_exception(
                  boost::format(
                    "unexpected character '%1%' in guide tree") % ch);
            break;
        }
        break;

      case st_Number:
        if (ch == '.')
          state = st_Fraction;
        else if (isalpha(ch))
          state = st_ID;
        else if (not isdigit(ch))
        {
          m_data.unget();
          m_token.erase(m_token.end() - 1);
          m_lookahead = gtt_Weight;
        }
        break;

      case st_Fraction:
        if (not isdigit(ch))
        {
          m_data.unget();
          m_token.erase(m_token.end() - 1);
          m_lookahead = gtt_Weight;
        }
        break;

      case st_ID:
        if (not isalnum(ch) and ch != '_')
        {
          m_data.unget();
          m_token.erase(m_token.end() - 1);
          m_lookahead = gtt_ID;
        }
        break;
    }
  }
}

void GuideTreeParser::match(GuideTreeToken token)
{
  if (token == m_lookahead)
    getNextToken();
  else
    throw mas_exception(
        boost::format(
          "invalid guide tree, expected %1% but found %2% ('%3%')") %
      int(token) % int(m_lookahead) % m_token);
}

tr1::tuple<base_node*,float> GuideTreeParser::parseNode()
{
  base_node* n = NULL;
  float w = 0;

  if (m_lookahead == gtt_Open)
    n = parseGroup();
  else
  {
    n = m_map[m_token];
    if (n == NULL)
      throw mas_exception(
          boost::format("guide tree contains unknown id %1%") % m_token);
    match(gtt_ID);
  }

  if (m_lookahead == gtt_Colon)
  {
    match(gtt_Colon);
    w = boost::lexical_cast<float>(m_token);
    match(gtt_Weight);
  }

  return tr1::make_tuple(n, w);
}

base_node* GuideTreeParser::parseGroup()
{
  base_node* na = NULL;
  base_node* nb = NULL;
  float wa = 0, wb = 0;

  match(gtt_Open);

  tr1::tie(na, wa) = parseNode();

  while (m_lookahead == gtt_Comma)
  {
    match(gtt_Comma);

    tr1::tie(nb, wb) = parseNode();

    na = new joined_node(na, nb, wa, wb);
    wa = wb;
  }

  match(gtt_Close);

  return na;
}

base_node* GuideTreeParser::parse()
{
  base_node* result = NULL;
  if (m_lookahead == gtt_Open)
    result = parseGroup();
  if (m_lookahead != gtt_End)
    throw mas_exception("invalid guide tree file, missing semicolon at end");
  return result;
}

void useGuideTree(const string& guide, vector<base_node*>& tree)
{
  map<string,leaf_node*> m;
  foreach (base_node* n, tree)
  {
    leaf_node* leaf = static_cast<leaf_node*>(n);
    m[leaf->m_entry.m_id] = leaf;
  }

  fs::ifstream file(guide);
  if (not file.is_open())
    throw mas_exception("failed to open guide tree");

  tree.clear();
  GuideTreeParser parser(file, m);
  tree.push_back(parser.parse());
}

