// Copyright Maarten L. Hekkelman, Radboud University 2008-2011.
// Copyright Coos Baakman, Jon Black, Wouter G. Touw & Gert Vriend, Radboud university medical center 2015.
//   Distributed under the Boost Software License, Version 1.0.
//       (See accompanying file LICENSE_1_0.txt or copy at
//             http://www.boost.org/LICENSE_1_0.txt)

#ifndef XSSP_ALIGN2D_H
#define XSSP_ALIGN2D_H

#pragma once

#include "mas.h"

#include <algorithm>
#include <string>
#include <vector>

typedef std::basic_string<uint8> sec_structure;

// --------------------------------------------------------------------
// entry is a multiple sequence alignment 'entry', a sequence with an ID and
// more.

struct entry
{
  entry(const entry& e)
    : m_nr(e.m_nr)
    , m_id(e.m_id)
    , m_seq(e.m_seq)
    , m_weight(e.m_weight)
    , m_positions(e.m_positions)
  {
  }

  entry(uint32 nr, const std::string& id, const sequence& seq = sequence(),
        float weight = 1.0f)
    : m_nr(nr)
    , m_id(id)
    , m_seq(seq)
    , m_weight(weight) {}

  int64 nr() const { return m_nr; }
  float weight() const { return m_weight; }
  uint32 length() const { return static_cast<uint32>(m_seq.length()); }

  void insert_gap(int64 pos);
  void append_gap();

  void remove_gap(int64 pos);

  void remove_gaps();
  void dump_positions() { m_positions.clear(); }

  int64 m_nr;
  std::string m_id;
  sequence m_seq;
  sec_structure m_ss;
  float m_weight;
  std::vector<int64> m_positions;
};

// --------------------------------------------------------------------
// nodes used in calculating a guide tree

struct base_node
{
  virtual        ~base_node() {}

  virtual void    print(std::ostream& s) = 0;

  virtual base_node*  left() const    { return 0; }
  virtual base_node*  right() const    { return 0; }

  virtual void    add_weight(float w) = 0;
  virtual uint32    leaf_count() const  { return 1; }

  virtual uint32    length() const = 0;
  virtual uint32    cost() const    { return 0; }
  virtual uint32    cumulative_cost() const
                      { return 0; }
};

std::ostream& operator<<(std::ostream& lhs, base_node& rhs);

struct joined_node : public base_node
{
  joined_node();
  joined_node(base_node* left, base_node* right, float d_left, float d_right);
  virtual ~joined_node();

  virtual void print(std::ostream& s);

  virtual base_node* left() const { return m_left; }
  virtual base_node* right() const { return m_right; }

  virtual void add_weight(float w)
  {
    m_left->add_weight(w);
    m_right->add_weight(w);
  }

  virtual uint32 leaf_count() const { return m_leaf_count; }
  virtual uint32 length() const { return m_length; }

  virtual uint32 cost() const { return m_length * m_leaf_count; }
  virtual uint32 cumulative_cost() const
  {
    return cost() + m_left->cumulative_cost() + m_right->cumulative_cost();
  }

  base_node* m_left;
  base_node* m_right;
  float m_d_left;
  float m_d_right;
  uint32 m_leaf_count;
  uint32 m_length;
};

struct leaf_node : public base_node
{
  leaf_node(entry& e)
    : m_entry(e)
  {
    m_entry.m_weight = 0;
  }

  virtual void print(std::ostream& s);

  virtual void add_weight(float w)
  {
    m_entry.m_weight += w;
  }

  virtual uint32 length() const
  {
    return static_cast<uint32>(m_entry.m_seq.length());
  }

  entry& m_entry;
};

class substitution_matrix_family;

template<typename T>
class matrix;

// prototype
void align(
  const joined_node* node,
  std::vector<entry*>& a, std::vector<entry*>& b, std::vector<entry*>& c,
  const substitution_matrix_family& mat_fam,
  float gop, float gep, float magic,
  bool ignorePositions);

void print_matrix(std::ostream& os,
  const matrix<int8>& tb, const sequence& sx, const sequence& sy);

#endif
