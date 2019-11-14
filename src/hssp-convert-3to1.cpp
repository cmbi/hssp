//  Copyright Maarten L. Hekkelman, Radboud University 2011.
// Copyright Coos Baakman, Jon Black, Wouter G. Touw & Gert Vriend, Radboud university medical center 2015.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "buffer.h"
#include "mas.h"
#include "utils.h"

#include <boost/algorithm/string.hpp>
#include <boost/filesystem/convenience.hpp>
#include <boost/filesystem/fstream.hpp>
#include <boost/foreach.hpp>
#include <boost/iostreams/copy.hpp>

#ifdef HAVE_LIBBZ2
#include <boost/iostreams/filter/bzip2.hpp>
#endif
#ifdef HAVE_LIBZ
#include <boost/iostreams/filter/gzip.hpp>
#endif

#include <boost/format.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/program_options.hpp>
#include <boost/program_options/config.hpp>
#include <boost/ptr_container/ptr_vector.hpp>
#include <boost/regex.hpp>
#include <boost/algorithm/string/regex.hpp>

#include <cmath>

#include "hssp-convert-3to1.h"

#if P_WIN
#pragma warning (disable: 4267)
#endif

#define foreach BOOST_FOREACH

namespace ba = boost::algorithm;
namespace io = boost::iostreams;
namespace fs = boost::filesystem;
namespace po = boost::program_options;

// --------------------------------------------------------------------
// basic named sequence type and a multiple sequence alignment container

struct insertion
{
  uint32 m_ipos, m_jpos;
  std::string m_seq;
};

class seq
{
  public:
    seq(const std::string& acc);
    seq(const seq&);
    ~seq();

  void validate(const seq& qseq);

  seq& operator=(const seq&);

  void swap(seq& o);

  std::string acc() const { return m_impl->m_acc; }

  void id(const std::string& id);
  std::string id() const { return m_impl->m_id; }

  void pdb(const std::string& pdb);
  std::string pdb() const { return m_impl->m_pdb; }

  void desc(const std::string& desc);
  std::string desc() const { return m_impl->m_desc; }

  void hssp(const std::string& hssp);

  float identity() const { return m_impl->m_identical; }
  float similarity() const { return m_impl->m_similar; }

  uint32 ifir() const { return m_impl->m_ifir; }
  uint32 ilas() const { return m_impl->m_ilas; }
  uint32 jfir() const { return m_impl->m_jfir; }
  uint32 jlas() const { return m_impl->m_jlas; }
  uint32 gapn() const { return m_impl->m_gapn; }
  uint32 gaps() const { return m_impl->m_gaps; }

  uint32 alignment_length() const { return m_impl->m_length; }
  uint32 seqlen() const { return m_impl->m_seqlen; }

  const std::list<insertion>& insertions() const
  {
    return m_impl->m_insertions;
  }

  void append(const std::string& seq);

  void update(const seq& qseq);
  static void  update_all(buffer<seq*>& b, const seq& qseq);

  bool operator<(const seq& o) const
  {
    return m_impl->m_identical > o.m_impl->m_identical or
           (m_impl->m_identical == o.m_impl->m_identical and
            length() > o.length());
  }

  uint32 length() const { return m_impl->m_end - m_impl->m_begin; }

  char& operator[](uint32 offset)
  {
    assert(offset < m_impl->m_size);
    return m_impl->m_seq[offset];
  }

  char operator[](uint32 offset) const
  {
    assert(offset < m_impl->m_size);
    return m_impl->m_seq[offset];
  }

  template<class T>
  class basic_iterator : public std::iterator<std::bidirectional_iterator_tag,T>
  {
  public:
    typedef typename std::iterator<std::bidirectional_iterator_tag, T> base_type;
    typedef typename base_type::reference reference;
    typedef typename base_type::pointer pointer;

    basic_iterator(T* s) : m_seq(s) {}
    basic_iterator(const basic_iterator& o) : m_seq(o.m_seq) {}

    basic_iterator& operator=(const basic_iterator& o)
    {
      m_seq = o.m_seq;
      return *this;
    }

    reference operator*() { return *m_seq; }
    reference operator->() { return *m_seq; }

    basic_iterator& operator++() { ++m_seq; return *this; }
    basic_iterator operator++(int)
    {
      basic_iterator iter(*this);
      operator++();
      return iter;
    }

    basic_iterator& operator--() { --m_seq; return *this; }
    basic_iterator operator--(int)
    {
      basic_iterator iter(*this);
      operator--();
      return iter;
    }

    bool operator==(const basic_iterator& o) const { return m_seq == o.m_seq; }
    bool operator!=(const basic_iterator& o) const { return m_seq != o.m_seq; }

    template<class U>
    friend basic_iterator<U> operator-(basic_iterator<U>, int);

    private:
      pointer m_seq;
  };

  typedef basic_iterator<char> iterator;
  typedef basic_iterator<const char> const_iterator;

  iterator begin() { return iterator(m_impl->m_seq); }
  iterator end() { return iterator(m_impl->m_seq + m_impl->m_size); }

  const_iterator begin() const
  {
    return const_iterator(m_impl->m_seq);
  }

  const_iterator end() const
  {
    return const_iterator(m_impl->m_seq + m_impl->m_size);
  }

  private:

    struct seq_impl
    {
      seq_impl(const std::string& acc);
      ~seq_impl();

      void update(const seq_impl& qseq);

      iterator begin() { return iterator(m_seq); }
      iterator end() { return iterator(m_seq + m_size); }

      const_iterator begin() const { return const_iterator(m_seq); }
      const_iterator end() const { return const_iterator(m_seq + m_size); }

      std::string m_id;
      std::string m_acc;
      std::string m_desc;
      std::string m_pdb;
      uint32 m_ifir;
      uint32 m_ilas;
      uint32 m_jfir;
      uint32 m_jlas;
      uint32 m_length;
      uint32 m_seqlen;
      float m_similar;
      float m_identical;
      uint32 m_begin, m_end;
      uint32 m_gaps, m_gapn;
      std::list<insertion> m_insertions;
      char* m_data;
      char* m_seq;
      uint32 m_refcount;
      uint32 m_size;
      uint32 m_space;
    };

    seq_impl* m_impl;

    seq();
};

template<class T>
inline seq::basic_iterator<T> operator-(seq::basic_iterator<T> i, int o)
{
  seq::basic_iterator<T> r(i);
  r.m_seq -= o;
  return r;
}

//typedef boost::ptr_std::vector<seq> mseq;
typedef std::vector<seq> mseq;

const uint32 kBlockSize = 512;

seq::seq_impl::seq_impl(const std::string& acc)
  : m_acc(acc),
    m_identical(0),
    m_similar(0),
    m_length(0),
    m_seqlen(0),
    m_begin(0),
    m_end(0),
    m_gaps(0),
    m_gapn(0),
    m_data(nullptr),
    m_seq(nullptr),
    m_refcount(1),
    m_size(0),
    m_space(0)
{
  m_ifir = m_ilas = m_jfir = m_jlas = 0;

  std::string::size_type s = m_acc.find('/');
  if (s != std::string::npos)
    m_acc.erase(s, std::string::npos);
}

seq::seq_impl::~seq_impl()
{
  assert(m_refcount == 0);
  delete m_data;
}

seq::seq(const std::string& acc)
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

void seq::id(const std::string& id)
{
  m_impl->m_id = id;
}

void seq::pdb(const std::string& pdb)
{
  m_impl->m_pdb = pdb;
}

void seq::desc(const std::string& desc)
{
  m_impl->m_desc = desc;
}

void seq::hssp(const std::string& hssp)
{
  // HSSP score=0.98/1.00 aligned=1-46/1-46 length=46 ngaps=0 gaplen=0
  // seqlen=46

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

void seq::append(const std::string& seq)
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

  //bool sgapped = false;
  bool qgapped = false;

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

    //if (sgap)
      //sgapped = true;
    //else if (qgap)
    if (qgap)
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
      if(not sgap) ++jpos;
    }
    else
    {
      if (qgapped)
      {
        *si = tolower(*si);
        ins.m_seq += *si;
        m_insertions.push_back(ins);
      }

      //sgapped = false;
      qgapped = false;

      ++ipos;
      if(not sgap) ++jpos;
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
  if (gaps != m_impl->m_gaps)
  {
    std::cerr << "gaps != m_gaps (" << gaps << " - " << m_impl->m_gaps << ')'
         << std::endl;
    error = true;
  }

  if (gapn != m_impl->m_gapn)\
  {
    std::cerr << "gapn != m_gapn (" << gapn << " - " << m_impl->m_gapn << ')'
         << std::endl;
    error = true;
  }

  if (len != m_impl->m_length)
  {
    std::cerr << "len != m_length (" << len << " - " << m_impl->m_length << ')'
         << std::endl;
    error = true;
  }

  if (m_impl->m_jfir + ylen - 1 != m_impl->m_jlas)
  {
    std::cerr << "jfir != jlas + jlen (" << m_impl->m_jfir << ", "
         << m_impl->m_jlas << ", " << ylen << ')' << std::endl;
    error = true;
  }

  // if (m_impl->m_ifir + xlen - 1 != m_impl->m_ilas)
  // {
  //   std::cerr << "ifir != ilas + ilen (" << m_impl->m_ifir << ", "
  //        << m_impl->m_ilas << ", " << xlen << ')' << std::endl;
  //   error = true;
  // }

  float score = boost::lexical_cast<float>(
      (boost::format("%4.2f") % (float(ident) / len)).str());

  if (abs(score - m_impl->m_identical) > 0.1)
  {
    std::cerr << "score != m_identical (" << score << ", "
         << m_impl->m_identical << ")" << std::endl;
    error = true;
  }

  if (error)
    throw mas_exception(
        boost::format("validation failed for %1%") % m_impl->m_id);
}

// --------------------------------------------------------------------

struct ResidueInfo
{
      ResidueInfo(const std::string& ri);

  std::string  m_ri, m_pr;
  uint32  m_pos;
};

ResidueInfo::ResidueInfo(const std::string& ri)
  : m_ri(ri), m_pos(0)
{
  for (int i = 34; i < 38; ++i)
    m_ri[i] = m_ri[i + 1];
  m_ri[38] = ' ';
}

typedef std::shared_ptr<ResidueInfo> res_ptr;
typedef std::vector<res_ptr>      res_list;

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

typedef std::shared_ptr<Hit> hit_ptr;
typedef std::vector<hit_ptr>  hit_list;

// --------------------------------------------------------------------

uint32 ReadHSSP2File(std::istream& is, std::string& id, std::string& header,
                     mseq& msa, hit_list& hits, res_list& residues,
                     uint32& nchain,
                     std::map<std::string, std::vector<std::string>> &sameChains)
{
  std::string line;
  std::string qid;

  uint32 offset = residues.size();
  uint32 result = 0;
  uint32 rix = offset;
  std::string::size_type ccOffset = 0, queryNr=0, idWidth = 0;
  char chainId = 0;
  boost::regex r1("Chain (.+) is considered to be the same as (.+(?:(?:, .+)* and .+)?)");
  std::map<std::string,uint32> index;

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
        std::string s = m[2].str();
        boost::sregex_token_iterator i(s.begin(), s.end(), boost::regex(" *(and|,) *"), -1), j;
        for (; i != j; i++)
        {
            sameChains[m[1]].push_back(*i);
            nchain ++;
        }

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
      uint32 nr = boost::lexical_cast<uint32>(
          ba::trim_copy(line.substr(8, 5))) - 1 + offset;

      if (nr >= residues.size())
        throw mas_exception("invalid input file");

      residues[nr]->m_pr = line.substr(14);
      continue;
    }

    if (ba::starts_with(line, "#=GS "))
    {
      line.erase(0, 5);
      // first GS line, fetch the width
      if (msa.size() == queryNr + 1 and ba::starts_with(line, qid))
      {
        ccOffset = 6 + line.find(" CC ");
      }
      else
      {
        std::string id = line.substr(0, ccOffset - 5);
        ba::trim(id);

        if (index.find(id) == index.end())
        {
          index[id] = msa.size();
          hits.push_back(
              hit_ptr(new Hit(msa.size(), queryNr, chainId, offset)));

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

      std::string id = line.substr(0, idWidth);
      ba::trim(id);

      if (id.length() <= 0)
        throw mas_exception(boost::format{
          "Encountered an empty id in alignment line: \"%1%\""} % line);

      std::string seq = line.substr(idWidth);

      if (id == qid)
      {
        if (index.find(id) == index.end())
          throw mas_exception(boost::format{
            "id \"%1%\" is encountered in the alignment, but not in the header"
          } % id);

        uint32 pos = msa[index[id]].length();
        foreach (char r, seq)
        {
          if (not is_gap(r))
          {
            ++result;

            if (rix < residues.size() and residues[rix]->m_ri[8] == '!')
              ++rix;

            if (r != residues[rix]->m_ri[8] and
                not (r == 'C' or islower(residues[rix]->m_ri[8])))
              throw mas_exception(boost::format {"Invalid hssp3 residue: \'%c\' expected to be \'%c\'"} % residues[rix]->m_ri[8] % r);

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
// Write collected information as a HSSP file to the output stream

void CreateHSSPOutput(const std::string& inProteinID,
                      const std::string& inProteinDescription,
                      float inThreshold, uint32 inSeqLength, uint32 inNChain,
                      uint32 inKChain, const std::map<std::string, std::vector<std::string>> &inSameChains,
                      const std::string& inUsedChains,
                      mseq& msa, hit_list& hits, res_list& res,
                      std::ostream& os)
{
  // print the header
  os << "HSSP       HOMOLOGY DERIVED SECONDARY STRUCTURE OF PROTEINS , VERSION " << PACKAGE_VERSION << " 2018" << std::endl
     << "PDBID      " << inProteinID << std::endl
     //<< "SEQBASE    " << inDatabank->GetName() << " version " << inDatabank->GetVersion() << std::endl
     << "THRESHOLD  according to: t(L)=(290.15 * L ** -0.562) + " << (inThreshold * 100) << std::endl
     << "REFERENCE  Sander C., Schneider R. : Database of homology-derived protein structures. Proteins, 9:56-68 (1991)." << std::endl
     << "CONTACT    Maintained at http://www.cmbi.umcn.nl/ <hssp.cmbi@radboudumc.nl>" << std::endl
     << inProteinDescription
     << boost::format("SEQLENGTH %5.5d") % inSeqLength << std::endl
     << boost::format("NCHAIN     %4.4d chain(s) in %s data set") % inNChain % inProteinID << std::endl;

  for (const auto &pair: inSameChains)
  {
    os << boost::format("SCHAIN     %s") % pair.first;
    for (const std::string &chain: pair.second)
      os << boost::format(", %s") % chain;
    os << std::endl;
  }

  if (inKChain != inNChain)
    os << boost::format("KCHAIN     %4.4d chain(s) used here ; chains(s) : ") % inKChain << inUsedChains << std::endl;

  os << boost::format("NALIGN     %4.4d") % hits.size() << std::endl
     << "NOTATION : ID: EMBL/SWISSPROT identifier of the aligned (homologous) protein" << std::endl
     << "NOTATION : STRID: if the 3-D structure of the aligned protein is known, then STRID is the Protein Data Bank identifier as taken" << std::endl
     << "NOTATION : from the database reference or DR-line of the EMBL/SWISSPROT entry" << std::endl
     << "NOTATION : %IDE: percentage of residue identity of the alignment" << std::endl
     << "NOTATION : %SIM (%WSIM):  (weighted) similarity of the alignment" << std::endl
     << "NOTATION : IFIR/ILAS: first and last residue of the alignment in the test sequence" << std::endl
     << "NOTATION : JFIR/JLAS: first and last residue of the alignment in the alignend protein" << std::endl
     << "NOTATION : LALI: length of the alignment excluding insertions and deletions" << std::endl
     << "NOTATION : NGAP: number of insertions and deletions in the alignment" << std::endl
     << "NOTATION : LGAP: total length of all insertions and deletions" << std::endl
     << "NOTATION : LSEQ2: length of the entire sequence of the aligned protein" << std::endl
     << "NOTATION : ACCNUM: SwissProt accession number" << std::endl
     << "NOTATION : PROTEIN: one-line description of aligned protein" << std::endl
     << "NOTATION : SeqNo,PDBNo,AA,STRUCTURE,BP1,BP2,ACC: sequential and PDB residue numbers, amino acid (lower case = Cys), secondary" << std::endl
     << "NOTATION : structure, bridge partners, solvent exposure as in DSSP (Kabsch and Sander, Biopolymers 22, 2577-2637(1983)" << std::endl
     << "NOTATION : VAR: sequence variability on a scale of 0-100 as derived from the NALIGN alignments" << std::endl
     << "NOTATION : pair of lower case characters (AvaK) in the alignend sequence bracket a point of insertion in this sequence" << std::endl
     << "NOTATION : dots (....) in the alignend sequence indicate points of deletion in this sequence" << std::endl
     << "NOTATION : SEQUENCE PROFILE: relative frequency of an amino acid type at each position. Asx and Glx are in their" << std::endl
     << "NOTATION : acid/amide form in proportion to their database frequencies" << std::endl
     << "NOTATION : NOCC: number of aligned sequences spanning this position (including the test sequence)" << std::endl
     << "NOTATION : NDEL: number of sequences with a deletion in the test protein at this position" << std::endl
     << "NOTATION : NINS: number of sequences with an insertion in the test protein at this position" << std::endl
     << "NOTATION : ENTROPY: entropy measure of sequence variability at this position" << std::endl
     << "NOTATION : RELENT: relative entropy, i.e.  entropy normalized to the range 0-100" << std::endl
     << "NOTATION : WEIGHT: conservation weight" << std::endl
     << std::endl
     << "## PROTEINS : identifier and alignment statistics" << std::endl
     << "  NR.    ID         STRID   %IDE %WSIM IFIR ILAS JFIR JLAS LALI NGAP LGAP LSEQ2 ACCNUM     PROTEIN" << std::endl;

  // print the first list
  uint32 nr = 1;
  boost::format fmt1("%5.5d : %12.12s%4.4s    %4.2f  %4.2f%5.5d%5.5d%5.5d%5.5d%5.5d%5.5d%5.5d%5.5d  %10.10s %s");
  foreach (hit_ptr h, hits)
  {
    const seq& s(msa[h->m_seq]);

    std::string id = s.id();
    if (id.length() > 12)
      id.erase(12, std::string::npos);
    else if (id.length() < 12)
      id.append(12 - id.length(), ' ');

    std::string acc = s.acc();
    if (acc.length() > 10)
      acc.erase(10, std::string::npos);
    else if (acc.length() < 10)
      acc.append(10 - acc.length(), ' ');

    std::string pdb = s.pdb();
    if (pdb.empty())
      pdb.append(4, ' ');

    os << fmt1 % nr
           % id % pdb
           % s.identity() % s.similarity()
           % (s.ifir() + h->m_offset) % (s.ilas() + h->m_offset)
           % s.jfir() % s.jlas() % s.alignment_length()
           % s.gaps() % s.gapn() % s.seqlen()
           % acc % s.desc()
       << std::endl;

    ++nr;
  }

  // print the alignments
  for (uint32 i = 0; i < hits.size(); i += 70)
  {
    uint32 m = i + 70,
           n = m,
           j;
    if (n > hits.size())
      n = hits.size();

    uint32 k[7] = {
      ((i +  0) / 10 + 1) % 10,
      ((i + 10) / 10 + 1) % 10,
      ((i + 20) / 10 + 1) % 10,
      ((i + 30) / 10 + 1) % 10,
      ((i + 40) / 10 + 1) % 10,
      ((i + 50) / 10 + 1) % 10,
      ((i + 60) / 10 + 1) % 10
    };

    os << boost::format("## ALIGNMENTS %4.4d - %4.4d") % (i + 1) % n << std::endl
       << boost::format(" SeqNo  PDBNo AA STRUCTURE BP1 BP2  ACC NOCC  VAR  ....:....%1.1d....:....%1.1d....:....%1.1d....:....%1.1d....:....%1.1d....:....%1.1d....:....%1.1d")
          % k[0] % k[1] % k[2] % k[3] % k[4] % k[5] % k[6] << " CHAIN AUTHCHAIN     NUMBER     RESNUM        BP1        BP2" << std::endl;

    res_ptr last;
    uint32 nr = 1;

    foreach (res_ptr ri, res)
    {
      std::string aln;

      if (ri->m_ri[6] != '!')
      {
        foreach (hit_ptr hit, boost::make_iterator_range(hits.begin() + i, hits.begin() + n))
        {
          const seq& s = msa[hit->m_seq];

          uint32 ifir = s.ifir() + hit->m_offset;
          uint32 ilas = s.ilas() + hit->m_offset;

          if (nr >= ifir and nr <= ilas)
            aln += msa[hit->m_seq][ri->m_pos];
          else
            aln += ' ';
        }

        // Add whitespace placeholders.
        for (j = n; j < m; j++)
          aln += ' ';
      }

      os << boost::format(" %5.5d") % nr << ri->m_ri.substr(0, 43) << "  " << aln << ri->m_ri.substr(43) << std::endl;
      ++nr;
    }
  }

  // ## SEQUENCE PROFILE AND ENTROPY
  os << "## SEQUENCE PROFILE AND ENTROPY" << std::endl
     << " SeqNo PDBNo   V   L   I   M   F   W   Y   G   A   P   S   T   C   H   R   K   Q   E   N   D  NOCC NDEL NINS ENTROPY RELENT WEIGHT CHAIN AUTHCHAIN     NUMBER     RESNUM        BP1        BP2" << std::endl;

  res_ptr last;
  nr = 1;
  foreach (res_ptr r, res)
  {
    if (r->m_pr.empty())
      os << boost::format("%5.5d") % nr << "          0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0     0    0    0   0.000      0  1.00" << std::endl;
    else
      os << boost::format("%5.5d") % nr << r->m_pr << std::endl;
    ++nr;
  }

  // insertion list

  os << "## INSERTION LIST" << std::endl
     << " AliNo  IPOS  JPOS   Len Sequence" << std::endl;

  foreach (hit_ptr h, hits)
  {
    seq& seq = msa[h->m_seq];

    //foreach (insertion& ins, h->insertions)
    foreach (const insertion& ins, seq.insertions())
    {
      std::string s = ins.m_seq;

      if (s.length() <= 100)
        os << boost::format(" %5.5d %5.5d %5.5d %5.5d ") % h->m_nr % (ins.m_ipos + h->m_offset) % ins.m_jpos % (ins.m_seq.length() - 2) << s << std::endl;
      else
      {
        os << boost::format(" %5.5d %5.5d %5.5d %5.5d ") % h->m_nr % (ins.m_ipos + h->m_offset) % ins.m_jpos % (ins.m_seq.length() - 2) << s.substr(0, 100) << std::endl;
        s.erase(0, 100);

        while (not s.empty())
        {
          uint32 n = s.length();
          if (n > 100)
            n = 100;

          os << "     +                   " << s.substr(0, n) << std::endl;
          s.erase(0, n);
        }
      }
    }
  }

  os << "//" << std::endl;
}

// --------------------------------------------------------------------

void ConvertHsspFile(std::istream& in, std::ostream& out)
{
  hit_list hits;
  res_list residues;

  std::string id, header;
  uint32 seqlength = 0, nchain = 0, kchain = 0;
  std::vector<std::string> usedChains;
  mseq msa;
  std::map<std::string, std::vector<std::string>> sameChains;


  for (;;)
  {
    std::string line;
    getline(in, line);
    if (line.empty() and in.eof())
      break;

    if (line != "# STOCKHOLM 1.0")
      throw mas_exception("Not a stockholm file, missing first line");

    if (not residues.empty())
      residues.push_back(res_ptr(new ResidueInfo("      ! !              0   0    0    0    0")));

    std::string h;
    swap(header, h);

    uint32 chainLength = ReadHSSP2File(in, id, header, msa, hits, residues,
                                       nchain, sameChains);
//    if (not h.empty() and h != header)
//      throw mas_exception("Inconsistent HSSP3 file, different header parts");

    if (chainLength == 0)
      break;

    ++kchain;
    seqlength += chainLength;
    usedChains.push_back(std::string(&hits.back()->m_chain, 1));
  }

  sort(hits.begin(), hits.end(),
    [&msa](hit_ptr a, hit_ptr b) -> bool
    {
      const seq& sa = msa[a->m_seq];
      const seq& sb = msa[b->m_seq];

      return sa.identity() > sb.identity() or
        (sa.identity() == sb.identity() and (
          (sa.id() < sb.id() or (sa.id() == sb.id()
                                 and sa.ifir() < sb.ifir()))));
    }
  );

  if (hits.size() > 9999)
    hits.erase(hits.begin() + 9999, hits.end());

  if (hits.empty())
    throw mas_exception("No hits found or remaining");

  uint32 nr = 1;
  foreach (hit_ptr h, hits)
    h->m_nr = nr++;

  CreateHSSPOutput(id, header, 0.05f, seqlength, nchain, kchain,
                   sameChains, ba::join(usedChains, ", "), msa,
                   hits, residues, out);
}
