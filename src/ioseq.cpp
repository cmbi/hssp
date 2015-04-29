// Copyright Maarten L. Hekkelman, Radboud University 2008-2011.
// Copyright Coos Baakman, Jon Black, Wouter G. Touw & Gert Vriend, Radboud university medical center 2015.
//   Distributed under the Boost Software License, Version 1.0.
//       (See accompanying file LICENSE_1_0.txt or copy at
//             http://www.boost.org/LICENSE_1_0.txt)
//
// i/o code for mas alignments and sequences

#include "ioseq.h"

#include "mas.h"
#include "structure.h"
#include "utils.h"

#include <boost/algorithm/string.hpp>
#include <boost/date_time/local_time/local_time.hpp>
#include <boost/filesystem/fstream.hpp>
#include <boost/foreach.hpp>
#include <boost/format.hpp>
#include <boost/iostreams/device/array.hpp>
#include <boost/iostreams/stream.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/regex.hpp>

#include <iostream>

using namespace std;
namespace fs = boost::filesystem;
namespace io = boost::iostreams;
namespace ba = boost::algorithm;

#define foreach BOOST_FOREACH
// --------------------------------------------------------------------

void readFasta(istream& is, vector<entry>& seq)
{
  string id, s;

  for (;;)
  {
    string line;
    getline(is, line);
    if (line.empty() or line[0] == '>')
    {
      if (line.empty() and not is.eof())
        continue;

      if (not (id.empty() or s.empty()))
      {
        if (VERBOSE)
          cout << "Sequence " << seq.size() + 1 << ": "
             << id << string(20 - id.length(), ' ') << s.length() << " aa" << endl;
        seq.push_back(entry(seq.size(), id, encode(s)));
      }
      id.clear();
      s.clear();

      if (line.empty())
        break;

      string::size_type w = line.find(' ');
      if (w != string::npos)
        id = line.substr(1, w - 1);
      else
        id = line.substr(1);

      s.clear();
      continue;
    }

    s += line;
  }
}

void readPDB(istream& is, char chainID, vector<entry>& seq)
{
  MProtein p(is);

  if (chainID == 0)
    chainID = p.GetFirstChainID();

  p.CalculateSecondaryStructure();

  entry e(seq.size(), p.GetID());
  p.GetSequence(chainID, e);
  seq.push_back(e);
}

void readAlignmentFromHsspFile(
  istream&    is,
  char&      chainID,
  vector<entry>&  seq)
{
  seq.clear();

  string line;

  // first line, should be something like "HSSP   bla bla bla"
  getline(is, line);
  if (not ba::starts_with(line, "HSSP"))
    throw mas_exception("file is not an HSSP file, does not start with HSSP");

  // second line, should contain PDBID
  getline(is, line);
  if (not ba::starts_with(line, "PDBID") or line.length() < 13)
    throw mas_exception("file is not a valid HSSP file, PDBID missing");

  entry pdb(0, line.substr(11), sequence());
  seq.push_back(pdb);

  string header;
  uint32 seqLength = 0, nchain = 0, nalign = 0, chain = 0;

  while (not is.eof())
  {
    getline(is, line);

    if (ba::starts_with(line, "## "))
      break;

    if (ba::starts_with(line, "HEADER     "))
      header = line.substr(11);
    else if (ba::starts_with(line, "SEQLENGTH  "))
      seqLength = boost::lexical_cast<uint32>(ba::trim_copy(line.substr(11, 4)));
    else if (ba::starts_with(line, "NCHAIN     "))
      nchain = boost::lexical_cast<uint32>(ba::trim_copy(line.substr(11, 4)));
    else if (ba::starts_with(line, "NALIGN     "))
      nalign = boost::lexical_cast<uint32>(ba::trim_copy(line.substr(11, 4)));
    else if (ba::starts_with(line, "KCHAIN     "))
    {
      vector<string> chains;
      string l = line.substr(48);
      ba::split(chains, l, ba::is_any_of(","));

      if (chains.empty())
        throw mas_exception("Invalid KCHAIN line in HSSP file");

      if (chainID == 0)
        chainID = chains.front()[0];
      else
      {
        for (chain = 0; chain < chains.size(); ++chain)
        {
          if (chainID == chains[chain][0])
            break;
        }

        if (chain == chains.size())
          throw mas_exception(boost::format("Chain %c not found in HSSP file") % chainID);
      }
    }
  }

  if (nalign < 1)
    throw mas_exception("invalid HSSP file, number of alignments missing");

  if (seqLength < 1)
    throw mas_exception("invalid HSSP file, sequence length missing");

  if (nchain < 1)
    throw mas_exception("invalid HSSP file, nchain missing");

  // second part, collect proteins

  if (line != "## PROTEINS : EMBL/SWISSPROT identifier and alignment statistics")
    throw mas_exception("invalid or unsupported HSSP file, ## PROTEINS line does match expected value");

  getline(is, line);
  if (not ba::starts_with(line, "  NR."))
    throw mas_exception("invalid HSSP file, expected line starting with NR.");

  seq.reserve(nalign);

  for (uint32 nr = 0; is.eof() == false and nr < nalign; ++nr)
  {
    getline(is, line);

    if (line.length() < 92)
      throw mas_exception("invalid HSSP file, protein line too short");

    entry e(nr + 1, line.substr(8, 12), sequence());
    ba::trim(e.m_id);

    string len = line.substr(74, 4);
    ba::trim(len);
    e.m_seq.reserve(boost::lexical_cast<uint32>(len));

    seq.push_back(e);
  }

  assert(seq.size() == nalign + 1);

  getline(is, line);

  uint32 alignment = 1;
  while (is.eof() == false and alignment < nalign)
  {
    if (not ba::starts_with(line, "## ALIGNMENTS"))
      throw mas_exception("invalid HSSP file, missing ## ALIGNMENTS line");

    string a_from = line.substr(14, 4);  ba::trim(a_from);
    string a_to = line.substr(21, 4); ba::trim(a_to);

    uint32 a_f = boost::lexical_cast<uint32>(a_from);
    uint32 a_t = boost::lexical_cast<uint32>(a_to);

    if (a_f != alignment or a_t < a_f or a_t > nalign)
      throw mas_exception("Invalid HSSP file, incorrect number of alignments");

    alignment = a_t + 1;

    getline(is, line);  // SeqNo line

    string pdb_seq;
    vector<string> s(a_t - a_f + 1);
    vector<int16> pdbnrs;

    do {
      getline(is, line);

      int32 n = static_cast<int32>(line.length()) - 51;

      if (chainID == 0)      // store chain ID
        chainID = line[12];

      if (line[12] == chainID)
      {
        string pdbno = line.substr(7, 4);
        ba::trim(pdbno);
        int16 pdbno_value = boost::lexical_cast<int16>(pdbno);

        pdbnrs.push_back(pdbno_value);

        pdb_seq += line[14];

        for (int32 i = 0; i < n; ++i)
        {
          char a = line[51 + i];
          if (a == ' ' or a == '.')
            a = '-';

          s[i] += a;
          seq[a_f + i].m_positions.push_back(pdbno_value);
        }
      }
    }
    while (not is.eof() and not (ba::starts_with(line, "##") or ba::starts_with(line, "//")));

    for (uint32 i = 0; i < a_t - a_f + 1; ++i)
    {
      seq[a_f + i].m_seq = encode(s[i]);
      assert(seq[a_f + i].m_seq.length() == seq[a_f + i].m_positions.size());
    }

    if (seq.front().m_seq.empty())
    {
      seq.front().m_seq = encode(pdb_seq);
      seq.front().m_positions = pdbnrs;
    }
    else if (decode(seq.front().m_seq) != pdb_seq)
    {
      cout << endl << decode(seq.front().m_seq) << endl
         << pdb_seq << endl
         << endl;

      throw mas_exception("Invalid HSSP file, inconsistent PDB sequence");
    }
  }

  // process ## INSERTION LIST, if any

  while (not is.eof() and ba::starts_with(line, "## ") and
    not ba::starts_with(line, "## INSERTION LIST"))
  {
    do
    {
      getline(is, line);
    }
    while (not (is.eof() or ba::starts_with(line, "## ")));
  }

  if (ba::starts_with(line, "## INSERTION LIST"))
  {
    getline(is, line);

    for (;;)
    {
      getline(is, line);

      if (is.eof() or ba::starts_with(line, "//"))
        break;

      uint32 l = boost::lexical_cast<uint32>(ba::trim_copy(line.substr(20, 4)));

      if (line.length() != 25 + l + 2)
        throw mas_exception("Invalid HSSP file, incorrect insertion length");

      sequence ins = encode(line.substr(26, l));

      uint32 p = boost::lexical_cast<uint32>(ba::trim_copy(line.substr(8, 4))) - 1;
      uint32 a = boost::lexical_cast<uint32>(ba::trim_copy(line.substr(2, 4)));

      seq[a].m_seq.insert(p, ins);
      seq[a].m_positions.insert(seq[a].m_positions.begin() + p, l, 0);
    }
  }

  seq.erase(
    remove_if(seq.begin(), seq.end(), boost::bind(&entry::length, _1) == 0UL),
    seq.end());

//  foreach (const entry& e, seq)
//  {
//    cout << e.m_id << endl
//       << decode(e.m_seq) << endl;
//
//    foreach (uint16 p, e.m_positions)
//      cout << p % 10;
//
//    cout << endl << endl;
//  }
}

void readWhatifMappingFile(istream& is, vector<entry>& seq)
{
  seq.clear();

  string line;

  while (not is.eof())
  {
    getline(is, line);
    if (not ba::starts_with(line, "Sequence name: "))
      continue;
    break;
  }

  while (ba::starts_with(line, "Sequence name: "))
  {
    string id = line.substr(15);
    string s;
    vector<int16> pos;

    // skip the description line
    getline(is, line);

    while (not (is.eof() or ba::starts_with(line, "Sequence name: ")))
    {
      getline(is, line);

      if (line.length() != 16 or line[4] != ' ' or line[6] != ' ' or line[11] != ' ')
        continue;

      s += line[5];

      string nr = line.substr(12);
      ba::trim(nr);

      if (nr == "----")
        pos.push_back(0);
      else
        pos.push_back(boost::lexical_cast<int16>(nr));
    }

    if (not s.empty())
    {
      entry e(seq.size(), id, encode(s));
      e.m_positions = pos;
      seq.push_back(e);
    }
  }
}

void readFamilyIdsFile(istream& is, const fs::path& dir, vector<entry>& seq)
{
  seq.clear();

  while (not is.eof())
  {
    string id;
    getline(is, id);

    if (id.empty())
      continue;

    fs::ifstream data(dir / (id + ".mapping"));
    if (not data.is_open())
      throw mas_exception(boost::format("Failed to open mapping file for protein %1%") % id);

    string line, s;
    vector<int16> pos;

    while (not data.eof())
    {
      getline(data, line);

      if (line.length() < 3 or line[1] != '\t')
        continue;

      s += line[0];
      pos.push_back(boost::lexical_cast<int16>(line.substr(2)));
    }

    if (not s.empty())
    {
      entry e(seq.size(), id, encode(s));
      e.m_positions = pos;
      seq.push_back(e);
    }
  }
}

void readSecStruct(std::vector<entry>& seq)
{
  throw mas_exception("please implement me first");
//  foreach (entry& e, seq)
//  {
//    fs::path ssfile(e.m_id + ".ss");
//    if (fs::exists(ssfile))
//    {
//      fs::ifstream ssdata(ssfile);
//
//      if (ssdata.is_open())
//      {
//        string line;
//        getline(ssdata, line);
//
//        if (encode(line) == e.m_seq)
//        {
//          getline(ssdata, line);
//          if (line.length() == e.m_seq.length())
//            e.m_ss = line;
//        }
//      }
//    }
//  }
}

// --------------------------------------------------------------------

void report_in_fasta(const vector<entry*>& alignment, ostream& os)
{
  foreach (const entry* e, alignment)
  {
    os << '>' << e->m_id << endl;

    uint32 o = 0;
    while (o < e->m_seq.length())
    {
      uint32 n = e->m_seq.length() - o;
      if (n > 72)
        n = 72;

      os << decode(e->m_seq.substr(o, n)) << endl;
      o += n;
    }
  }
}

void report_in_clustalw(const vector<entry*>& alignment, ostream& os)
{
  os << "CLUSTAL FORMAT for MaartensAlignment" << endl;

  uint32 len = alignment[0]->m_seq.length();
  uint32 offset = 0;

//  if (alignment.size() == 2)
//  {
//    // first strip off leading and trailing unaligned seqs
//
//    if (alignment.front()->m_seq[0] == kSignalGapCode)
//    {
//      do ++offset;
//      while (offset < len and alignment.front()->m_seq[offset] == kSignalGapCode);
//    }
//    else if (alignment.back()->m_seq[0] == kSignalGapCode)
//    {
//      do ++offset;
//      while (offset < len and alignment.back()->m_seq[offset] == kSignalGapCode);
//    }
//
//    if (*(alignment.front()->m_seq.begin() + len - 1) == kSignalGapCode)
//    {
//      do --len;
//      while (len > offset and *(alignment.front()->m_seq.begin() + len - 1) == kSignalGapCode);
//    }
//    else if (*(alignment.back()->m_seq.begin() + len - 1) == kSignalGapCode)
//    {
//      do --len;
//      while (len > offset and *(alignment.back()->m_seq.begin() + len - 1) == kSignalGapCode);
//    }
//  }

  while (offset < len)
  {
    uint32 n = len - offset;
    if (n > 60)
      n = 60;

    struct {
      uint32    cnt[20];
    } dist[60] = {};

    foreach (const entry* e, alignment)
    {
      string id = e->m_id;
      if (id.length() > 15)
        id = id.substr(0, 12) + "...";
      else if (id.length() < 15)
        id += string(15 - id.length(), ' ');

      if (VERBOSE > 1 and not e->m_ss.empty())
      {
        sec_structure s2 = e->m_ss.substr(offset, n);
        const char kSS[] = " HBEGITS";
        os << id << ' ';
        for (sec_structure::iterator i = s2.begin(); i != s2.end(); ++i)
          os << kSS[*i];
        os << endl;
      }

      sequence ss = e->m_seq.substr(offset, n);

      for (uint32 i = 0; i < n; ++i)
      {
        aa ri = ss[i];
        if (ri < 20)
          dist[i].cnt[ri] += 1;
      }

      os << id << ' ' << decode(ss) << endl;
    }

    string scores(n, ' ');
    for (uint32 i = 0; i < n; ++i)
    {
      uint32 strong[9] = {};
      const char* kStrongGroups[9] = {
        "STA", "NEQK", "NHQK", "NDEQ", "QHRK", "MILV", "MILF", "HY", "FYW"
      };

      uint32 weak[11] = {};
      const char* kWeakGroups[11] = {
        "CSA", "ATV", "SAG", "STNK", "STPA", "SGND", "SNDEQK",
        "NDEQHK", "NEQHRK", "FVLIM", "HFY"
      };

      for (uint32 r = 0; r < 20; ++r)
      {
        if (dist[i].cnt[r] == alignment.size())
        {
          scores[i] = '*';
          break;
        }

        for (uint32 g = 0; g < 9; ++g)
        {
          if (strchr(kStrongGroups[g], kAA[r]) != NULL)
            strong[g] += dist[i].cnt[r];
        }

        for (uint32 g = 0; g < 11; ++g)
        {
          if (strchr(kWeakGroups[g], kAA[r]) != NULL)
            weak[g] += dist[i].cnt[r];
        }
      }

      for (uint32 g = 0; scores[i] == ' ' and g < 9; ++g)
      {
        if (strong[g] == alignment.size())
          scores[i] = ':';
      }

      for (uint32 g = 0; scores[i] == ' ' and g < 11; ++g)
      {
        if (weak[g] == alignment.size())
          scores[i] = '.';
      }
    }

    os << string(16, ' ') << scores << endl;

    if (not alignment.front()->m_positions.empty() and VERBOSE > 1)
    {
      string pos_nrs(n, ' ');
      for (uint32 i = 0; i < n; ++i)
      {
        if (alignment.front()->m_positions[offset + i] != 0)
          pos_nrs[i] = '!';
      }

      os << string(16, ' ') << pos_nrs << endl;
    }

    offset += n;
    os << endl;
  }
}

uint32 gcg_checksum(const sequence& seq)
{
  uint32 result = 0, i = 0;

  foreach (aa r, seq)
    result += ((i++ % 57) + 1) * kAA[r];

  return result % 10000;
}

void report_in_msf(const vector<entry*>& alignment, ostream& os)
{
  using namespace boost::posix_time;

  uint32 len = alignment.front()->m_seq.length();

  os << "!!AA_MULTIPLE_ALIGNMENT 1.0" << endl
     << endl
     << " Mas MSF: " << len
     << " Type: P " << second_clock::local_time()
     << " Check: 0" << endl
     << endl;

  foreach (const entry* e, alignment)
  {
    os << "Name: " << e->m_id << ' '
       << "Len: " << e->m_seq.length() << ' '
       << "Check: " << gcg_checksum(e->m_seq) << ' '
       << "Weight: " << e->m_weight
       << endl;
  }

  os << endl
     << "//" << endl
     << endl;

  uint32 offset = 0;
  while (offset < len)
  {
    uint32 n = len - offset;
    if (n > 50)
      n = 50;

    string l1 = boost::lexical_cast<string>(offset + 1);
    string l2 = boost::lexical_cast<string>(offset + 50);

    os << string(22, ' ') << l1 << string(54 - l1.length() - l2.length(), ' ') << l2 << endl;

    foreach (const entry* e, alignment)
    {
      string ss = decode(e->m_seq.substr(offset, n));
      for (uint32 i = 4; i > 0; --i)
      {
        if (ss.length() > i * 10)
          ss.insert(ss.begin() + i * 10, ' ');
      }

      string id = e->m_id;
      if (id.length() > 21)
        id = id.substr(0, 18) + "...";
      else if (id.length() < 21)
        id += string(21 - id.length(), ' ');

      os << id << ' ' << ss << endl;
    }

    offset += n;
    os << endl;
  }
}

void report(const vector<entry*>& alignment, ostream& os, const string& format)
{
  if (format == "fasta")
    report_in_fasta(alignment, os);
  else if (format == "clustalw")
    report_in_clustalw(alignment, os);
  else if (format == "msf")
    report_in_msf(alignment, os);
  else
    throw mas_exception(boost::format("Unknown output format %1%") % format);
}

