//  Copyright Maarten L. Hekkelman, Radboud University 2011.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include "buffer.h"
#include "mas.h"
#include "matrix.h"
#include "utils.h"

#include <boost/algorithm/string.hpp>
#include <boost/filesystem/convenience.hpp>
#include <boost/filesystem/fstream.hpp>
#include <boost/foreach.hpp>
#include <boost/format.hpp>
#include <boost/iostreams/copy.hpp>
#include <boost/iostreams/filter/bzip2.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/program_options.hpp>
#include <boost/program_options/config.hpp>
#include <boost/ptr_container/ptr_vector.hpp>
#include <boost/regex.hpp>

#include <cmath>
#include <iostream>
#include <numeric>

#if P_WIN
#pragma warning (disable: 4267)
#endif

using namespace std;
namespace ba = boost::algorithm;
namespace io = boost::iostreams;
namespace fs = boost::filesystem;
namespace po = boost::program_options;

#define foreach BOOST_FOREACH
int VERBOSE = 0;

// --------------------------------------------------------------------
// basic named sequence type and a multiple sequence alignment container

struct insertion
{
  uint32      m_ipos, m_jpos;
  string      m_seq;
};

class seq
{
  public:
        seq(const string& acc);
        seq(const seq&);
        ~seq();

  void    validate(const seq& qseq);

  seq&    operator=(const seq&);

  void    swap(seq& o);

  string    acc() const              { return m_impl->m_acc; }

  void    id(const string& id);
  string    id() const              { return m_impl->m_id; }

  void    pdb(const string& pdb);
  string    pdb() const              { return m_impl->m_pdb; }

  void    desc(const string& desc);
  string    desc() const            { return m_impl->m_desc; }

  void    hssp(const string& hssp);

  float    identity() const          { return m_impl->m_identical; }
  float    similarity() const          { return m_impl->m_similar; }

  uint32    ifir() const            { return m_impl->m_ifir; }
  uint32    ilas() const            { return m_impl->m_ilas; }
  uint32    jfir() const            { return m_impl->m_jfir; }
  uint32    jlas() const            { return m_impl->m_jlas; }
  uint32    gapn() const            { return m_impl->m_gapn; }
  uint32    gaps() const            { return m_impl->m_gaps; }

  uint32    alignment_length() const      { return m_impl->m_length; }
  uint32    seqlen() const            { return m_impl->m_seqlen; }

  const list<insertion>&
        insertions() const          { return m_impl->m_insertions; }

  void    append(const string& seq);

  void    update(const seq& qseq);
  static void  update_all(buffer<seq*>& b, const seq& qseq);

  bool    operator<(const seq& o) const    { return m_impl->m_identical > o.m_impl->m_identical or
                            (m_impl->m_identical == o.m_impl->m_identical and length() > o.length()); }

  uint32    length() const            { return m_impl->m_end - m_impl->m_begin; }

  char&    operator[](uint32 offset)
        {
          assert(offset < m_impl->m_size);
          return m_impl->m_seq[offset];
        }

  char    operator[](uint32 offset) const
        {
          assert(offset < m_impl->m_size);
          return m_impl->m_seq[offset];
        }

  template<class T>
  class basic_iterator : public std::iterator<bidirectional_iterator_tag,T>
  {
    public:
    typedef typename std::iterator<std::bidirectional_iterator_tag, T>  base_type;
    typedef  typename base_type::reference                reference;
    typedef typename base_type::pointer                  pointer;

            basic_iterator(T* s) : m_seq(s) {}
            basic_iterator(const basic_iterator& o) : m_seq(o.m_seq) {}

    basic_iterator&  operator=(const basic_iterator& o)
            {
              m_seq = o.m_seq;
              return *this;
            }

    reference    operator*()          { return *m_seq; }
    reference    operator->()        { return *m_seq; }

    basic_iterator&  operator++()        { ++m_seq; return *this; }
    basic_iterator  operator++(int)        { basic_iterator iter(*this); operator++(); return iter; }

    basic_iterator&  operator--()        { --m_seq; return *this; }
    basic_iterator  operator--(int)        { basic_iterator iter(*this); operator--(); return iter; }

    bool      operator==(const basic_iterator& o) const
                          { return m_seq == o.m_seq; }
    bool      operator!=(const basic_iterator& o) const
                          { return m_seq != o.m_seq; }

    template<class U>
    friend basic_iterator<U> operator-(basic_iterator<U>, int);

    private:
    pointer      m_seq;
  };

  typedef basic_iterator<char>    iterator;
  typedef basic_iterator<const char>  const_iterator;

  iterator    begin()              { return iterator(m_impl->m_seq); }
  iterator    end()              { return iterator(m_impl->m_seq + m_impl->m_size); }

  const_iterator  begin() const          { return const_iterator(m_impl->m_seq); }
  const_iterator  end() const            { return const_iterator(m_impl->m_seq + m_impl->m_size); }

  private:

  struct seq_impl
  {
          seq_impl(const string& acc);
          ~seq_impl();

    void    update(const seq_impl& qseq);

    iterator  begin()              { return iterator(m_seq); }
    iterator  end()              { return iterator(m_seq + m_size); }

    const_iterator
          begin() const          { return const_iterator(m_seq); }
    const_iterator
          end() const            { return const_iterator(m_seq + m_size); }

    string    m_id, m_acc, m_desc, m_pdb;
    uint32    m_ifir, m_ilas, m_jfir, m_jlas, m_length, m_seqlen;
    float    m_identical, m_similar;
    uint32    m_begin, m_end;
    uint32    m_gaps, m_gapn;
    list<insertion>
          m_insertions;
    char*    m_data;
    char*    m_seq;
    uint32    m_refcount;
    uint32    m_size, m_space;
  };

  seq_impl*  m_impl;

        seq();
};

template<class T>
inline seq::basic_iterator<T> operator-(seq::basic_iterator<T> i, int o)
{
  seq::basic_iterator<T> r(i);
  r.m_seq -= o;
  return r;
}

//typedef boost::ptr_vector<seq> mseq;
typedef vector<seq>        mseq;

const uint32 kBlockSize = 512;

seq::seq_impl::seq_impl(const string& acc)
  : m_acc(acc)
  , m_length(0)
  , m_identical(0)
  , m_similar(0)
  , m_seqlen(0)
  , m_begin(0)
  , m_end(0)
  , m_gaps(0)
  , m_gapn(0)
  , m_data(nullptr)
  , m_seq(nullptr)
  , m_refcount(1)
  , m_size(0)
  , m_space(0)
{
  m_ifir = m_ilas = m_jfir = m_jlas = 0;

  string::size_type s = m_acc.find('/');
  if (s != string::npos)
    m_acc.erase(s, string::npos);
}

seq::seq_impl::~seq_impl()
{
  assert(m_refcount == 0);
  delete m_data;
}

seq::seq(const string& acc)
  : m_impl(new seq_impl(acc))
{
}

seq::seq(const seq& s)
  : m_impl(s.m_impl)
{
  ++m_impl->m_refcount;
}

seq& seq::operator=(const seq& rhs)
{
  if (this != &rhs)
  {
    if (--m_impl->m_refcount == 0)
      delete m_impl;

    m_impl = rhs.m_impl;

    ++m_impl->m_refcount;
  }

  return *this;
}

seq::~seq()
{
  if (--m_impl->m_refcount == 0)
    delete m_impl;
}

void seq::swap(seq& o)
{
  std::swap(m_impl, o.m_impl);
}

void seq::id(const string& id)
{
  m_impl->m_id = id;
}

void seq::pdb(const string& pdb)
{
  m_impl->m_pdb = pdb;
}

void seq::desc(const string& desc)
{
  m_impl->m_desc = desc;
}

void seq::hssp(const string& hssp)
{
  // HSSP score=0.98/1.00 aligned=1-46/1-46 length=46 ngaps=0 gaplen=0 seqlen=46

  static const boost::regex
    re1("score=(\\d\\.\\d+)/(\\d\\.\\d+)"),
    re2("aligned=(\\d+)-(\\d+)/(\\d+)-(\\d+)"),
    re3("length=(\\d+)"),
    re4("ngaps=(\\d+)"),
    re5("gaplen=(\\d+)"),
    re6("seqlen=(\\d+)");

  boost::smatch m;
  if (boost::regex_search(hssp, m, re1))
  {
    m_impl->m_identical = boost::lexical_cast<float>(m[1]);
    m_impl->m_similar = boost::lexical_cast<float>(m[2]);
  }

  if (boost::regex_search(hssp, m, re2))
  {
    m_impl->m_ifir = boost::lexical_cast<uint32>(m[1]);
    m_impl->m_ilas = boost::lexical_cast<uint32>(m[2]);
    m_impl->m_jfir = boost::lexical_cast<uint32>(m[3]);
    m_impl->m_jlas = boost::lexical_cast<uint32>(m[4]);
  }

  if (boost::regex_search(hssp, m, re3))
    m_impl->m_length = boost::lexical_cast<uint32>(m[1]);

  if (boost::regex_search(hssp, m, re4))
    m_impl->m_gaps = boost::lexical_cast<uint32>(m[1]);

  if (boost::regex_search(hssp, m, re5))
    m_impl->m_gapn = boost::lexical_cast<uint32>(m[1]);

  if (boost::regex_search(hssp, m, re6))
    m_impl->m_seqlen = boost::lexical_cast<uint32>(m[1]);
}

void seq::append(const string& seq)
{
  if (m_impl->m_size + seq.length() > m_impl->m_space)
  {
    // increase storage for the sequences
    uint32 k = m_impl->m_space;
    if (k == 0)
      k = kBlockSize;
    uint32 n = k * 2;
    if (n < seq.length())
      n = seq.length();
    char* p = new char[n];
    memcpy(p, m_impl->m_data, m_impl->m_size);
    delete [] m_impl->m_data;
    m_impl->m_data = m_impl->m_seq = p;
    m_impl->m_space = n;
  }

  memcpy(m_impl->m_seq + m_impl->m_size, seq.c_str(), seq.length());
  m_impl->m_end = m_impl->m_size += seq.length();
}

void seq::update_all(buffer<seq*>& b, const seq& qseq)
{
  for (;;)
  {
    seq* s = b.get();
    if (s == nullptr)
      break;

    s->update(qseq);
  }

  b.put(nullptr);
}

void seq::update(const seq& qseq)
{
  m_impl->update(*qseq.m_impl);
}

void seq::seq_impl::update(const seq_impl& qseq)
{
  uint32 ipos = 1, jpos = m_jfir;
  if (jpos == 0)
    jpos = 1;

  bool sgapped = false, qgapped = false;

  const_iterator qi = qseq.begin();
  iterator si = begin();
  uint32 i = 0;
  insertion ins = {};

  for (; qi != qseq.end(); ++qi, ++si, ++i)
  {
    bool qgap = is_gap(*qi);
    bool sgap = is_gap(*si);

    if (qgap and sgap)
      continue;

    if (sgap)
      sgapped = true;
    else if (qgap)
    {
      if (not qgapped)
      {
        iterator gsi = si - 1;
        while (gsi != begin() and is_gap(*gsi))
          --gsi;

        ins.m_ipos = ipos;
        ins.m_jpos = jpos;
        ins.m_seq = *gsi = tolower(*gsi);
      }

      ins.m_seq += *si;

      qgapped = true;
      ++jpos;
    }
    else
    {
      if (qgapped)
      {
        *si = tolower(*si);
        ins.m_seq += *si;
        m_insertions.push_back(ins);
      }

      sgapped = false;
      qgapped = false;

      ++ipos;
      ++jpos;
    }
  }
}

void seq::validate(const seq& qseq)
{
  uint32 gaps = 0, gapn = 0, len = 0, xlen = 0, ylen = 0, ident = 0;

  if (qseq.length() != length())
    throw mas_exception("Invalid sequence length");

  bool started = false;
  bool xgapped = false;
  bool ygapped = false;

  const char* xseq = qseq.m_impl->m_seq;
  const char* yseq = m_impl->m_seq;

  uint32 n = m_impl->m_size;
  while (n > 0 and (is_gap(xseq[n - 1]) or is_gap(yseq[n - 1])))
    --n;

  for (uint32 i = 0; i < n; ++i)
  {
    bool xgap = is_gap(xseq[i]);
    bool ygap = is_gap(yseq[i]);

    if (not started and (xgap or ygap))
      continue;

    started = true;

    if (xgap and ygap)
      continue;

    ++len;
    if (xseq[i] == toupper(yseq[i]))
      ++ident;

    if (xgap)
    {
      if (not xgapped)
        ++gaps;
      ++gapn;
      xgapped = true;
    }
    else
    {
      ++xlen;
      xgapped = false;
    }

    if (ygap)
    {
      if (not ygapped)
        ++gaps;
      ++gapn;
      ygapped = true;
    }
    else
    {
      ++ylen;
      ygapped = false;
    }
  }

  bool error = false;
  if (gaps != m_impl->m_gaps)      { cerr << "gaps != m_gaps (" << gaps << " - " << m_impl->m_gaps << ')' << endl; error = true; }
  if (gapn != m_impl->m_gapn)      { cerr << "gapn != m_gapn (" << gapn << " - " << m_impl->m_gapn << ')' << endl; error = true; }
  if (len != m_impl->m_length)    { cerr << "len != m_length (" << len << " - " << m_impl->m_length << ')' << endl; error = true; }
  if (m_impl->m_jfir + ylen - 1 != m_impl->m_jlas)
                    { cerr << "jfir != jlas + jlen (" << m_impl->m_jfir << ", " << m_impl->m_jlas << ", " << ylen << ')' << endl; error = true; }
//  if (m_impl->m_ifir + xlen - 1 != m_impl->m_ilas)
//                    { cerr << "ifir != ilas + ilen (" << m_impl->m_ifir << ", " << m_impl->m_ilas << ", " << xlen << ')' << endl; error = true; }

  float score = boost::lexical_cast<float>((boost::format("%4.2f") % (float(ident) / len)).str());

  if (abs(score - m_impl->m_identical) > 0.1) { cerr << "score != m_identical (" << score << ", " << m_impl->m_identical << ")" << endl; error = true; }

  if (error)
    throw mas_exception(boost::format("validation failed for %1%") % m_impl->m_id);
}

namespace std
{
  template<>
  void swap(seq& a, seq& b)
  {
    a.swap(b);
  }
}

// --------------------------------------------------------------------

struct ResidueInfo
{
      ResidueInfo(const string& ri);

  string  m_ri, m_pr;
  uint32  m_pos;
};

ResidueInfo::ResidueInfo(const string& ri)
  : m_ri(ri), m_pos(0)
{
  for (int i = 34; i < 38; ++i)
    m_ri[i] = m_ri[i + 1];
  m_ri[38] = ' ';
}

typedef shared_ptr<ResidueInfo> res_ptr;
typedef vector<res_ptr>      res_list;

// --------------------------------------------------------------------

struct Hit
{
          Hit(uint32 s, uint32 q, char chain, uint32 offset)
            : m_seq(s), m_qseq(q), m_chain(chain), m_nr(0), m_offset(offset)
          {
          }

  uint32      m_seq, m_qseq;
  char      m_chain;
  uint32      m_nr, m_offset;
};

typedef shared_ptr<Hit> hit_ptr;
typedef vector<hit_ptr>  hit_list;

// --------------------------------------------------------------------

uint32 ReadHSSP2File(istream& is, string& id, string& header, mseq& msa, hit_list& hits, res_list& residues, uint32& nchain)
{
  string line, qid;

  uint32 offset = residues.size(), result = 0, rix = offset;
  string::size_type ccOffset = 0, queryNr, idWidth = 0;
  char chainId;
  boost::regex r1("Chain . is considered to be the same as (.(?:(?:, .)* and .)?)");
  map<string,uint32> index;

  nchain += 1;

  for (;;)
  {
    line.clear();
    getline(is, line);

    if (line.empty())
    {
      if (not is.good())
        throw mas_exception("Stockholm file is truncated or incomplete");
      continue;
    }

    if (line == "//")
      break;

    if (ba::starts_with(line, "#=GF ID "))
    {
      queryNr = msa.size();
      qid = line.substr(8);
      index[qid] = queryNr;
      msa.push_back(seq(qid));
      continue;
    }

    if (ba::starts_with(line, "#=GF CC "))
    {
      line.erase(0, 8);

      if (ba::starts_with(line, "PDBID "))
      {
        id = line.substr(7);
        continue;
      }

      if (ba::starts_with(line, "DATE   "))
      {
        header += line.substr(0, 7) + "    file generated on " + line.substr(7) + '\n';
        continue;
      }

      if (ba::starts_with(line, "PDBID  ") or
        ba::starts_with(line, "HEADER ") or
        ba::starts_with(line, "COMPND ") or
        ba::starts_with(line, "SOURCE ") or
        ba::starts_with(line, "AUTHOR ") or
        ba::starts_with(line, "DBREF  "))
      {
        header += line.substr(0, 7) + "    " + line.substr(7) + '\n';
        continue;
      }

      boost::smatch m;

      if (boost::regex_match(line, m, r1))
      {
        string s = m[1];
        if (s.length() == 1)
          nchain += 1;
        else
          nchain += (s.length() - 1) / 3;

        continue;
      }

      continue;
    }

    if (ba::starts_with(line, "#=GF RI "))
    {
      chainId = line[20];
      residues.push_back(res_ptr(new ResidueInfo(line.substr(14))));
      continue;
    }

    if (ba::starts_with(line, "#=GF PR "))
    {
      uint32 nr = boost::lexical_cast<uint32>(ba::trim_copy(line.substr(8, 5))) - 1 + offset;
      if (nr >= residues.size())
        throw mas_exception("invalid input file");

      residues[nr]->m_pr = line.substr(14);
      continue;
    }

    if (ba::starts_with(line, "#=GS "))
    {
      line.erase(0, 5);
      if (msa.size() == queryNr + 1 and ba::starts_with(line, qid))  // first GS line, fetch the width
      {
        ccOffset = 6 + line.find(" CC ");
      }
      else
      {
        string id = line.substr(0, ccOffset - 5);
        ba::trim(id);

        if (index.find(id) == index.end())
        {
          index[id] = msa.size();
          hits.push_back(hit_ptr(new Hit(msa.size(), queryNr, chainId, offset)));
          msa.push_back(seq(id));
        }

        line.erase(0, ccOffset - 5);

        if (ba::starts_with(line, "ID "))
          msa[index[id]].id(line.substr(3));
        else if (ba::starts_with(line, "DE "))
          msa[index[id]].desc(line.substr(3));
        else if (ba::starts_with(line, "HSSP "))
          msa[index[id]].hssp(line.substr(5));
        else if (ba::starts_with(line, "DR PDB "))
          msa[index[id]].pdb(line.substr(7, 4));
      }

      continue;
    }

    if (line[0] != '#')
    {
      if (idWidth == 0)
      {
        assert(ba::starts_with(line, qid));
        idWidth = qid.length();
        while (idWidth < line.length() and line[idWidth] == ' ')
          ++idWidth;
      }

      string id = line.substr(0, idWidth);
      ba::trim(id);

      string seq = line.substr(idWidth);

      if (id == qid)
      {
        uint32 pos = msa[index[id]].length();
        foreach (char r, seq)
        {
          if (not is_gap(r))
          {
            ++result;

            if (rix < residues.size() and residues[rix]->m_ri[8] == '!')
              ++rix;

            if (r != residues[rix]->m_ri[8] and not (r == 'C' or islower(residues[rix]->m_ri[8])))
              throw mas_exception("Invalid hssp3 file");

            residues[rix]->m_pos = pos;
            ++rix;
          }

          ++pos;
        }
      }

      msa[index[id]].append(seq);
    }
  }

  for (uint32 i = queryNr + 1; i < msa.size(); ++i)
    msa[i].update(msa[queryNr]);

//  for (uint32 i = queryNr + 1; i < msa.size(); ++i)
//    msa[i].validate(msa[queryNr]);

  return result;
}

// --------------------------------------------------------------------

class Alignment
{
  public:
        Alignment(const mseq& msa);
        ~Alignment();

  uint32    M() const                { return m_m; }
  uint32    N() const                { return m_n; }

  uint8&    operator()(uint32 i, uint32 j)      { assert(i < m_m); assert(j < m_n); return m_data[i * m_n + j]; }
  uint8    operator()(uint32 i, uint32 j) const  { assert(i < m_m); assert(j < m_n); return m_data[i * m_n + j]; }

  private:

  const mseq&  m_msa;
  uint8*    m_data;
  uint32    m_m, m_n;
};

Alignment::Alignment(const mseq& msa)
  : m_msa(msa)
  , m_data(nullptr)
{
  m_m = msa.size();
  m_n = 0;

  for (uint32 i = 0; i < msa[0].length(); ++i)
  {
    if (is_gap(msa[0][i]))
      continue;
    ++m_n;
  }

  m_data = new uint8[m_n * m_m];

  for (uint32 i = 0, j = 0; i < msa[0].length(); ++i)
  {
    if (is_gap(msa[0][i]))
      continue;

    for (uint32 k = 0; k < m_m; ++k)
    {
      uint8 r = ResidueNr(msa[k][i]);
      if (r > 20)
        r = 20;
      operator()(k, j) = r;
    }
    ++j;
  }
}

Alignment::~Alignment()
{
  delete[] m_data;
}

inline uint32 mapkey(uint32 i, uint32 alpha, uint32 q)
{
  return (q - 1) * (i - 1) + alpha;
}

void CalculateMI(const mseq& msa, ostream& os)
{
  Alignment align(msa);

  const double theta = 0.2;
  const double pseudocount_weight = 0.5;

  uint32 M = align.M();
  uint32 N = align.N();
  const uint32 q = 21;

  vector<double> W(M, 1.0f);

  if (theta > 0)
  {
    symmetric_matrix<bool> pdist(M);
    for (uint32 i = 0; i + 1 < M; ++i)
    {
      for (uint32 j = i + 1; j < M; ++j)
      {
        uint32 d = 0;
        for (uint32 k = 0; k < N; ++k)
        {
          if (align(i, k) != align(j, k))
            ++d;
        }

        pdist(i, j) = (static_cast<double>(d) / N) < theta;
      }
    }

    for (uint32 i = 0; i < M; ++i)
    {
      uint32 w = 1;
      for (uint32 j = 0; j < M; ++j)
        w += pdist(i, j);

      W[i] = 1.0f / w;
    }
  }

  double Meff = accumulate(W.begin(), W.end(), 0.f);

  matrix<double> Pi(N, q, 0.f);

  for (uint32 j = 0; j < M; ++j)
  {
    for (uint32 i = 0; i < N; ++i)
      Pi(i, align(j, i)) += W[j];
  }

  Pi /= Meff;

  const matrix<double> Q0(q, q);
  symmetric_matrix<matrix<double>> Pij(N, Q0);

  for (uint32 l = 0; l < M; ++l)
  {
    for (uint32 i = 0; i + 1 < N; ++i)
    {
      for (uint32 j = i + 1; j < N; ++j)
      {
        Pij(i, j)(align(l, i), align(l, j)) += W[l];
      }
    }
  }

  Pij /= Meff;
/*
  for (uint32 i = 0; i < N; ++i)
  {
    for (uint32 alpha = 0; alpha < q; ++alpha)
    {
      for (uint32 beta = 0; beta < q; ++beta)
      {
        if (alpha == beta)
          Pij(i, i)(alpha, beta) = Pi(i, alpha);
        else
          Pij(i, i)(alpha, beta) = 0;
      }
    }
  }
*/
  symmetric_matrix<double> MI(N);

  for (uint32 i = 0; i + 1 < N; ++i)
  {
    for (uint32 j = i + 1; j < N; ++j)
    {
      double m = 0;

      for (uint32 alpha = 0; alpha < q; ++alpha)
      {
        for (uint32 beta = 0; beta < q; ++beta)
        {
          if (Pij(i, j)(alpha, beta) > 0)
            m += Pij(i, j)(alpha, beta) * log(Pij(i, j)(alpha, beta) / (Pi(i, alpha) * Pi(j, beta)));
        }
      }

      MI(i, j) = m;
    }
  }

  os << "N: " << align.N() << endl
     << "M: " << align.M() << endl
     << "Meff: " << Meff << endl
     << endl;

  os << "  ";
  for (uint32 i = 0; i < N; ++i)
    os << kResidues[align(0, i)] << "     ";
  os << endl;

  for (uint32 i = 0; i < N; ++i)
  {
    os << kResidues[align(0, i)];

    for (uint32 j = 0; j < N; ++j)
      os << boost::format(" %5.3f") % MI(i, j);

    os << endl;
  }

  os << "//" << endl;

  // direct information

  // use pseudo count weight

  Pij *= 1 - pseudocount_weight;
  Pij += pseudocount_weight / (q * q);

  Pi *= 1 - pseudocount_weight;
  Pi += pseudocount_weight / q;

  for (uint32 i = 0; i < N; ++i)
  {
    for (uint32 alpha = 0; alpha < q; ++alpha)
    {
      for (uint32 beta = 0; beta < q; ++beta)
      {
        Pij(i, i)(alpha, beta) =
          Pij(i, i)(alpha, beta) - (pseudocount_weight / (q * q)) +
            (alpha == beta ? pseudocount_weight / q : 0);
      }
    }
  }

  // calculate inverse Covariance matrix

  symmetric_matrix<double> C(N);
  for (uint32 i = 0; i + 1 < N; ++i)
  {
    for (uint32 j = i + 1; j < N; ++j)
    {
      for (uint32 alpha = 0; alpha < q; ++alpha)
      {
        for (uint32 beta = 0; beta < q; ++beta)
        {
          C(mapkey(i, alpha, q), mapkey(j, beta, q)) =
            Pij(i, j)(alpha, beta) - Pi(i, alpha) * Pi(j, beta);
        }
      }
    }
  }

  C = invert(C);

  symmetric_matrix<double> DI(N);

  for (uint32 i = 0; i + 1 < N; ++i)
  {
    for (uint32 j = i + 1; j < N; ++j)
    {
      // direct information from mean-field

      matrix<double> W_mf(q, q, 1.0);

      for (uint32 ki = 0; ki < q - 1; ++ki)
      {
        for (uint32 kj = 0; kj < q - 1; ++kj)
        {
          W_mf(ki, kj) = exp(-C(mapkey(i, ki, q), mapkey(j, kj, q)));
        }
      }

      // calculate mu

      double mu1 = 1.0 / q, mu2 = 1.0 / q;
      double epsilon = 1e-4;
      double diff = 1.0;



      for (uint32 alpha = 0; alpha < q; ++alpha)
      {
        for (uint32 beta = 0; beta < q; ++beta)
        {
          if (Pij(i, j)(alpha, beta) > 0)
            m += Pij(i, j)(alpha, beta) * log(Pij(i, j)(alpha, beta) / (Pi(i, alpha) * Pi(j, beta)));
        }
      }

      MI(i, j) = m;
    }
  }



//  os << "  ";
//  for (uint32 i = 0; i < N; ++i)
//    os << kResidues[align(0, i)] << "     ";
//  os << endl;
//
//  for (uint32 i = 0; i < N; ++i)
//  {
//    os << kResidues[align(0, i)];
//
//    for (uint32 j = 0; j < N; ++j)
//      os << boost::format(" %5.3f") % MI(i, j);
//
//    os << endl;
//  }
//
//  os << "//" << endl;


}

// --------------------------------------------------------------------

void CalculateMI(istream& in, ostream& out)
{
  hit_list hits;
  res_list residues;

  string id, header;
  uint32 seqlength = 0, nchain = 0, kchain = 0;
  mseq msa;

  for (;;)
  {
    string line;
    getline(in, line);
    if (line.empty() and in.eof())
      break;

    if (line != "# STOCKHOLM 1.0")
      throw mas_exception("Not a stockholm file, missing first line");

    if (not residues.empty())
      residues.push_back(res_ptr(new ResidueInfo("      ! !              0   0    0    0    0")));
    uint32 first = residues.size();

    string h;
    swap(header, h);

    uint32 chainLength = ReadHSSP2File(in, id, header, msa, hits, residues, nchain);
    if (h.empty())
      out << header << endl;
    else if (h != header)
      throw mas_exception("Inconsistent HSSP3 file, different header parts");

    out << "Mutual information for chain " << residues.back()->m_ri[6] << endl;

    CalculateMI(msa, out);

    if (chainLength == 0)
      break;

    ++kchain;
    seqlength += chainLength;
  }
}

// --------------------------------------------------------------------

int main(int argc, char* const argv[])
{
  try
  {
    po::options_description desc("MKHSSP options");
    desc.add_options()
      ("help,h",               "Display help message")
      ("input,i",    po::value<string>(), "Input HSSP file")
      ("output,o",  po::value<string>(), "Output file, use 'stdout' to output to screen")
      ;

    po::positional_options_description p;
    p.add("input", 1);
    p.add("output", 2);

    po::variables_map vm;
    po::store(po::command_line_parser(argc, argv).options(desc).positional(p).run(), vm);

    //fs::path home = get_home();
    //if (fs::exists(home / ".mkhssprc"))
    //{
    //  fs::ifstream rc(home / ".mkhssprc");
    //  po::store(po::parse_config_file(rc, desc), vm);
    //}

    po::notify(vm);

    if (vm.count("help"))
    {
      cerr << desc << endl;
      exit(1);
    }

    io::filtering_stream<io::input> in;
    fs::ifstream ifs;

    if (vm.count("input") == 0)
      in.push(cin);
    else
    {
      fs::path input = vm["input"].as<string>();
      ifs.open(input, ios::binary);

      if (not ifs.is_open())
        throw mas_exception(boost::format("Could not open input file '%s'") % input);

      if (input.extension() == ".bz2")
        in.push(io::bzip2_decompressor());
      else if (input.extension() == ".gz")
        in.push(io::gzip_decompressor());
      in.push(ifs);
    }

    if (vm.count("output") == 0)
      CalculateMI(in, cout);
    else
    {
      fs::path output = vm["output"].as<string>();

      fs::ofstream ofs(output, ios::binary);
      if (not ofs.is_open())
        throw mas_exception(boost::format("Could not open output file '%s'") % output);

      io::filtering_stream<io::output> out;

      if (output.extension() == ".bz2")
        out.push(io::bzip2_compressor());
      else if (output.extension() == ".gz")
        out.push(io::gzip_compressor());
      out.push(ofs);

      CalculateMI(in, out);
    }
  }
  catch (const exception& e)
  {
    cerr << e.what() << endl;
    exit(1);
  }
  catch (...)
  {
    cerr << "Unknown exception" << endl;
    exit(1);
  }

  return 0;
}
