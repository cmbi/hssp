#include "hssp-nt.h"

#include "blast.h"
#include "buffer.h"
#ifdef HAVE_LIBZEEP
  #include "fetchdbrefs.h"
#endif
#include "matrix.h"
#include "progress.h"
#include "structure.h"
#include "utils.h"

#include <boost/algorithm/string.hpp>
#include <boost/date_time/date_clock_device.hpp>
#include <boost/date_time/gregorian/gregorian.hpp>
#include <boost/filesystem/convenience.hpp>
#include <boost/filesystem/fstream.hpp>
#include <boost/filesystem/operations.hpp>
#include <boost/foreach.hpp>
#include <boost/format.hpp>
#include <boost/iostreams/copy.hpp>
#include <boost/iostreams/device/back_inserter.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/range/adaptor/sliced.hpp>
#include <boost/regex.hpp>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <set>

#define foreach BOOST_FOREACH

namespace fs = boost::filesystem;
namespace ba = boost::algorithm;
namespace io = boost::iostreams;

namespace HSSP
{

// --------------------------------------------------------------------

const float kThreshold = 0.05f, kFragmentCutOff = 0.75f;

// precalculated threshold table for identity values between 10 and 80
const float kHomologyThreshold[] = {
  0.795468f, 0.75398f, 0.717997f, 0.686414f, 0.658413f, 0.633373f, 0.610811f,
  0.590351f, 0.571688f, 0.554579f, 0.53882f, 0.524246f, 0.510718f, 0.498117f,
  0.486344f, 0.475314f, 0.464951f, 0.455194f, 0.445984f, 0.437275f, 0.429023f,
  0.421189f, 0.413741f, 0.406647f, 0.399882f, 0.39342f, 0.38724f, 0.381323f,
  0.375651f, 0.370207f, 0.364976f, 0.359947f, 0.355105f, 0.35044f, 0.345941f,
  0.341599f, 0.337406f, 0.333352f, 0.329431f, 0.325636f, 0.32196f, 0.318396f,
  0.314941f, 0.311587f, 0.308331f, 0.305168f, 0.302093f, 0.299103f, 0.296194f,
  0.293362f, 0.290604f, 0.287917f, 0.285298f, 0.282744f, 0.280252f, 0.277821f,
  0.275448f, 0.273129f, 0.270865f, 0.268652f, 0.266488f, 0.264372f, 0.262302f,
  0.260277f, 0.258294f, 0.256353f, 0.254452f, 0.252589f, 0.250764f, 0.248975f,
  0.247221f
};

// --------------------------------------------------------------------

float calculateDistance(const sequence& a, const sequence& b)
{
  const float kDistanceGapOpen = 10;
  const float kDistanceGapExtend = 0.2f;

  int32 x = 0, dimX = static_cast<int32>(a.length());
  int32 y = 0, dimY = static_cast<int32>(b.length());

  matrix<float>  B(dimX, dimY);
  matrix<float>  Ix(dimX, dimY);
  matrix<float>  Iy(dimX, dimY);
  matrix<uint16>  id(dimX, dimY);

  Ix(0, 0) = 0;
  Iy(0, 0) = 0;

  uint16 highId = 0;

  Ix(x, y) = 0;
  Iy(x, y) = 0;
  if (x > 0 and y > 0)
    id(x - 1, y - 1) = highId;

  int32 startX = x, startY = y;
  float high = -std::numeric_limits<float>::max();
  uint16 highIdSub = 0;

  sequence::const_iterator ia = a.begin();
  for (x = startX; x < dimX; ++x, ++ia)
  {
    sequence::const_iterator ib = b.begin();
    for (y = startY; y < dimY; ++y, ++ib)
    {
      float Ix1 = 0; if (x > startX) Ix1 = Ix(x - 1, y);
      float Iy1 = 0; if (y > startY) Iy1 = Iy(x, y - 1);

      // (1)
      float M = score(kMPam250, *ia, *ib);
      if (x > startX and y > startY)
        M += B(x - 1, y - 1);

      float s;
      uint32 i = 0;
      if (a[x] == b[y])
        i = 1;

      if (M >= Ix1 and M >= Iy1)
      {
        if (x > startX and y > startY)
          i += id(x - 1, y - 1);
        s = M;
      }
      else if (Ix1 >= Iy1)
      {
        if (x > startX)
          i += id(x - 1, y);
        s = Ix1;
      }
      else
      {
        if (y > startY)
          i += id(x, y - 1);
        s = Iy1;
      }

      B(x, y) = s;
      id(x, y) = i;

      if ((x == dimX - 1 or y == dimY - 1) and high < s)
      {
        high = s;
        highIdSub = i;
      }

      // (3)
      Ix(x, y) = std::max(M - kDistanceGapOpen, Ix1 - kDistanceGapExtend);

      // (4)
      Iy(x, y) = std::max(M - kDistanceGapOpen, Iy1 - kDistanceGapExtend);
    }
  }

  highId += highIdSub;

  float result = 1.0f - float(highId) / std::max(dimX, dimY);

  assert(result >= 0.0f);
  assert(result <= 1.0f);

  return result;
}

// ----------------------------------------------------

char GetAACode(char oneLetterCode, int64 bridgeNumberSS)
{
    if (bridgeNumberSS != 0)
        return 'a' + (bridgeNumberSS - 1) % 26;
    else
        return oneLetterCode;
}

// ----------------------------------------------------

std::string GetStructureString(const MResidue &residue)
{
    std::string s(9, ' ');

    switch (residue.GetSecondaryStructure())
    {
        case alphahelix:  s[0] = 'H'; break;
        case betabridge:  s[0] = 'B'; break;
        case strand:      s[0] = 'E'; break;
        case helix_3:     s[0] = 'G'; break;
        case helix_5:     s[0] = 'I'; break;
        case turn:        s[0] = 'T'; break;
        case bend:        s[0] = 'S'; break;
        case loop:        s[0] = ' '; break;
    }

    for (uint32 stride = 3; stride <= 5; ++stride)
    {
        switch (residue.GetHelixFlag(stride))
        {
            case helixNone:         s[stride - 1] = ' '; break;
            case helixStart:        s[stride - 1] = '>'; break;
            case helixEnd:          s[stride - 1] = '<'; break;
            case helixStartAndEnd:  s[stride - 1] = 'X'; break;
            case helixMiddle:       s[stride - 1] = '0' + stride; break;
        }
    }

    double alpha;
    char chirality;
    std::tie(alpha, chirality) = residue.Alpha();

    s[5] = residue.IsBend()? 'S' : ' ';
    s[6] = chirality;

    for (uint32 i = 0; i < 2; ++i)
    {
        MBridgeParner p = residue.GetBetaPartner(i);
        if (p.residue != nullptr)
        {
            s[7 + i] = 'A' + p.ladder % 26;
            if (p.parallel)
                s[7 + i] = tolower(s[7 + i]);
        }
    }

    return s;
}


// --------------------------------------------------------------------

struct MResInfo
{
  uint8 m_letter;
  std::string m_chain_id,
              m_auth_chain_id;
  int64 m_seq_nr;
  int64 m_pdb_nr;
  std::string m_ins_code;
  MSecondaryStructure m_ss;
  char m_structure[10];
  char m_aa;
  int64 m_beta_partner_1,
        m_beta_partner_2,
        m_ss_bridge_nr;
  double m_accessibility;
  float m_consweight;
  uint32 m_nocc;
  uint32 m_dist[23];
  float m_dist_weight[23];
  float m_sum_dist_weight;
  uint32 m_ins, m_del;
  float m_score[23];
  float m_freq[20];
  float m_entropy;

  void Add(uint8 r, float inDistance);
  static MResInfo NewGap(size_t inDel, float inSumDistance, uint8 inResidue,
                         float inDistance);
  void AddGap(float inDistance);
};

typedef std::vector<MResInfo> MResInfoList;

void MResInfo::Add(uint8 r, float inDistance)
{
  float weight = 1 - inDistance;

  assert(r < 23);
  if (r < 22)
    m_nocc += 1;
  m_dist[r] += 1;

  m_dist_weight[r] += weight;
  m_sum_dist_weight += weight;

  for (int i = 0; i < 23; ++i)
  {
    float si = 0;

    for (int j = 0; j < 23; ++j)
      si += score(kMPam250, i, j) * m_dist_weight[j];

    m_score[i] = si / m_sum_dist_weight;
  }
}

void MResInfo::AddGap(float inDistance)
{
  Add(22, inDistance);
}

MResInfo MResInfo::NewGap(size_t inDel, float inSumDistance, uint8 inResidue,
                          float inDistance)
{
  MResInfo r = {};

  r.m_ins_code = " ";
  r.m_dist[22] = static_cast<uint32>(inDel);
  r.m_sum_dist_weight = r.m_dist_weight[22] = inSumDistance;
  r.Add(inResidue, inDistance);

  return r;
}

// --------------------------------------------------------------------

struct MHit;
typedef std::shared_ptr<MHit> MHitPtr;

struct MHit
{
  MHit(const MHit& e);

  MHit(const std::string& id, const std::string& def, const sequence& seq)
    : m_id(id),
      m_def(def),
      m_seq(seq),
      m_distance(0),
      m_identical(0),
      m_similar(0),
      m_length(0),
      m_gaps(0),
      m_gapn(0),
      m_score(0.0),
      m_ifir(0),
      m_ilas(0),
      m_jfir(0),
      m_jlas(0)
  {
  }

  struct insertion
  {
    uint32 m_ipos;
    uint32 m_jpos;
    std::string m_seq;
  };

  static MHitPtr Create(const std::string& id, const std::string& def,
                        const std::string& seq);

  void CalculateDistance(const sequence& chain);

  std::string m_id, m_acc, m_def, m_stid;
  sequence m_seq;
  std::string m_aligned;
  float m_distance, m_score;
  int32 m_ifir, m_ilas, m_jfir, m_jlas;
  uint32 m_identical, m_similar, m_length;
  uint32 m_gaps, m_gapn;
  std::vector<insertion> m_insertions;
};

std::ostream& operator<<(std::ostream& os, const MHit& hit)
{
  std::string seq = hit.m_aligned;
  foreach (char& r, seq)
    if (r == ' ' or r == '.') r = '-';

  for (std::string::size_type i = 72; i < seq.length(); i += 73)
    seq.insert(seq.begin() + i, '\n');

  os << '>' << hit.m_id /*<< ' ' << hit.m_def*/ << std::endl
     << seq << std::endl;

  return os;
}

MHitPtr MHit::Create(const std::string& id, const std::string& def,
                     const std::string& seq)
{
  MHitPtr result(new MHit(id, def, encode(seq)));

  static const boost::regex
    kM6FastARE("^(\\w+)((?:\\|([^| ]*))(?:\\|([^| ]+))?(?:\\|([^| ]+))?(?:\\|([^| ]+))?)");

  boost::smatch m;
  if (boost::regex_match(result->m_id, m, kM6FastARE,
                         boost::match_not_dot_newline))
  {
    if (m[1] == "sp" or m[1] == "tr")
    {
      result->m_acc = m[3];
      result->m_id = m[4];
    }
    else
      result->m_id = result->m_acc = m[2];
  }
  else
    result->m_acc = result->m_id;

  result->m_distance = 0;

  return result;
}

void MHit::CalculateDistance(const sequence& chain)
{
  m_distance = calculateDistance(m_seq, chain);
}

// --------------------------------------------------------------------

struct MProfile
{
  MProfile(const MChain& inChain, const sequence& inSequence,
           float inThreshold, float inFragmentCutOff);
  ~MProfile();

  void Process(std::istream& inHits, float inGapOpen, float inGapExtend,
               uint32 inMaxHits, uint32 inThreads);
  void Align(MHitPtr e, float inGapOpen, float inGapExtend);

  void AdjustXGapCosts(std::vector<float>& gop, std::vector<float>& gep);
  void AdjustYGapCosts(const sequence& s, std::vector<float>& gop,
                       std::vector<float>& gep);
  void dump(const matrix<float>& B, const matrix<float>& Ix,
            const matrix<float>& Iy, const matrix<int8>& tb,
            const std::vector<float>& gopX, const std::vector<float>& gopY,
            const std::vector<float>& gepX, const std::vector<float>& gepY,
            const sequence& sx, const sequence& sy);

  void PrintStockholm(std::ostream& os, const std::string& inChainID,
                      bool inFetchDBRefs) const;
  void PrintStockholm(std::ostream& os, const MProtein& inProtein,
                      bool inFetchDBRefs,
                      const std::vector<std::string>& inUsed,
                      const std::vector<std::string>& inAKA) const;

  void CalculateConservation(uint32 inThreads);

  const MChain&  m_chain;
  sequence m_seq;
  MResInfoList m_residues;
  std::vector<MHitPtr> m_entries;
  float m_threshold;
  float m_frag_cutoff;
  float m_sum_dist_weight;
  bool m_shuffled;
};

MProfile::MProfile(const MChain& inChain, const sequence& inSequence,
                   float inThreshold, float inFragmentCutOff)
  : m_chain(inChain),
    m_seq(inSequence),
    m_threshold(inThreshold),
    m_frag_cutoff(inFragmentCutOff),
    m_sum_dist_weight(0),
    m_shuffled(false)
{
  const std::vector<MResidue*>& residues = m_chain.GetResidues();
  std::vector<MResidue*>::const_iterator ri = residues.begin();

  int64 seq_nr = 1;
  for (uint32 i = 0; i < inSequence.length(); ++i)
  {
    assert(ri != residues.end());

    if (ri != residues.begin() and
        (*ri)->GetNumber() > (*(ri - 1))->GetNumber() + 1)
      ++seq_nr;

    int64 bridge_nr = 0;
    if ((*ri)->GetType() == kCysteine)
        bridge_nr = (*ri)->GetSSBridgeNr();

    MResInfo res = {
        inSequence[i], m_chain.GetChainID(), m_chain.GetAuthChainID(), seq_nr,
        (*ri)->GetSeqNumber(), (*ri)->GetInsertionCode(), (*ri)->GetSecondaryStructure()
    };
    res.m_aa = kResidueInfo[(*ri)->GetType()].code;
    strncpy(res.m_structure, GetStructureString(*(*ri)).c_str(), 9);
    res.m_beta_partner_1 = ((*ri)->GetBetaPartner(0).residue != nullptr)? (*ri)->GetBetaPartner(0).residue->GetNumber() : 0;
    res.m_beta_partner_2 = ((*ri)->GetBetaPartner(1).residue != nullptr)? (*ri)->GetBetaPartner(1).residue->GetNumber() : 0;
    res.m_ss_bridge_nr = bridge_nr;
    res.m_accessibility = (*ri)->Accessibility();
    res.Add(res.m_letter, 0);
    m_residues.push_back(res);

    ++ri;
    ++seq_nr;
  }
}

MProfile::~MProfile()
{
}

void MProfile::dump(const matrix<float>& B, const matrix<float>& Ix,
                    const matrix<float>& Iy, const matrix<int8>& tb,
                    const std::vector<float>& gopX,
                    const std::vector<float>& gopY,
                    const std::vector<float>& gepX,
                    const std::vector<float>& gepY,
                    const sequence& sx, const sequence& sy)
{
  std::ofstream os("alignment.log");
  os.imbue(std::locale(""));

  assert(sx.length() == m_residues.size());
  assert(sx.length() == gopX.size());
  assert(sx.length() == B.dim_m());
  assert(gopY.size() == B.dim_n());

  for (uint32 x = 0; x < sx.length(); ++x)
    os << '\t' << (is_gap(sx[x]) ? '.' : kResidues[sx[x]]) << '\t' << gopX[x];
  os << std::endl;

  for (uint32 y = 0; y < sy.length(); ++y)
  {
    os << kResidues[sy[y]];
    for (uint32 x = 0; x < m_residues.size(); ++x)
      os << '\t' << B(x, y) << '\t' << Iy(x, y);
    os << std::endl
       << gopY[y];

    for (uint32 x = 0; x < m_residues.size(); ++x)
    {
      switch (tb(x, y))
      {
        case -1:  os << '\t' << Ix(x, y) << '\t' << "|"; break;
        case  0:  os << '\t' << Ix(x, y) << '\t' << "\\"; break;
        case  1:  os << '\t' << Ix(x, y) << '\t' << "-"; break;
        case  2:  os << '\t' << Ix(x, y) << '\t' << "."; break;
      }
    }
    os << std::endl;
  }
}

void MProfile::AdjustXGapCosts(std::vector<float>& gop,
                               std::vector<float>& gep)
{
  assert(gop.size() == m_seq.length());
  assert(gop.size() == m_residues.size());

  for (size_t ix = 0; ix < m_residues.size(); ++ix)
  {
    MResInfo& e = m_residues[ix];

    // adjust for secondary structure
    switch (e.m_ss)
    {
      case alphahelix:  gop[ix] *= 2; break;
      case betabridge:  gop[ix] *= 2; break;
      case strand:    break;
      case helix_3:    gop[ix] *= 3; break;
      case helix_5:    gop[ix] *= 2; break;
      case turn:      gop[ix] *= 1.5; break;
      case bend:      break;
      case loop:      break;
    }

    // if there is a gap in the alignments, lower gap penalties
    // gap open penalty is zero when another gap was already created here
    if (e.m_ins > 0)
      gop[ix] = 0;

    // lower gap extension penalty for existing gaps here
    if (e.m_dist[22] > 0)
      gep[ix] = float(e.m_dist[22]) / (m_entries.size() + 1);

    // if there is a gap within 8 residues, increase gap penalty
    for (size_t d = 0; d < 8; ++d)
    {
      if (ix + d >= m_residues.size() or
        m_residues[ix + d].m_dist[22] > 0 or
        ix < d or
        m_residues[ix - d].m_dist[22] > 0)
      {
        gop[ix] *= (2 + ((8 - d) * 2)) / 8.f;
        break;
      }
    }
  }
}

const float kResidueSpecificPenalty[22] = {
  1.13f,    // A
  1.13f,    // C
  0.96f,    // D
  1.31f,    // E
  1.20f,    // F
  0.61f,    // G
  1.00f,    // H
  1.32f,    // I
  0.96f,    // K
  1.21f,    // L
  1.29f,    // M
  0.63f,    // N
  0.74f,    // P
  1.07f,    // Q
  0.72f,    // R
  0.76f,    // S
  0.89f,    // T
  1.25f,    // V
  1.23f,    // W
  1.00f,    // Y
  1.00f,    // B
  1.00f,    // Z
};

void MProfile::AdjustYGapCosts(const sequence& s, std::vector<float>& gop,
                               std::vector<float>& gep)
{
  for (uint32 y = 0; y < s.length(); ++y)
  {
    if (s[y] < 22)
      gop[y] *= kResidueSpecificPenalty[s[y]];
  }
}

void MProfile::Align(MHitPtr e, float inGapOpen, float inGapExtend)
{
  int32 x = 0, dimX = static_cast<int32>(m_seq.length());
  int32 y = 0, dimY = static_cast<int32>(e->m_seq.length());

  matrix<float> B(dimX, dimY);
  matrix<float> Ix(dimX, dimY);
  matrix<float> Iy(dimX, dimY);
  matrix<int8> tb(dimX, dimY);

  float minLength = static_cast<float>(dimX);
  float maxLength = static_cast<float>(dimY);
  if (minLength > maxLength)
    std::swap(minLength, maxLength);

  float logmin = 1.0f / log10(minLength);
  float logdiff = 1.0f + 0.5f * log10(minLength / maxLength);

  // initial gap open penalty
  float gop = (inGapOpen / (logdiff * logmin)) * abs(kMPam250MisMatchAverage) * kMPam250ScalingFactor;
  float gep = inGapExtend;

  // position specific gap penalties
  // initial gap extend penalty is adjusted for difference in sequence lengths
  std::vector<float> gop_a(dimX, gop);
  std::vector<float> gep_a(dimX, gep * (1 + log10(float(dimX) / dimY)));
  AdjustXGapCosts(gop_a, gep_a);

  std::vector<float> gop_b(dimY, gop);
  std::vector<float> gep_b(dimY, gep * (1 + log10(float(dimY) / dimX)));
  AdjustYGapCosts(e->m_seq, gop_b, gep_b);

  int32 highX = 0, highY = 0;
  float highS = 0;

  for (x = 0; x < dimX; ++x)
  {
    for (y = 0; y < dimY; ++y)
    {
      float Ix1 = 0; if (x > 0) Ix1 = Ix(x - 1, y);
      float Iy1 = 0; if (y > 0) Iy1 = Iy(x, y - 1);

      float M = m_residues[x].m_score[e->m_seq[y]];
      if (x > 0 and y > 0)
        M += B(x - 1, y - 1);

      if (M >= Ix1 and M >= Iy1)
      {
        tb(x, y) = 0;
        B(x, y) = M;

        Ix(x, y) = M - (x < dimX - 1 ? gop_a[x] : 0);
        Iy(x, y) = M - (y < dimY - 1 ? gop_b[y] : 0);

        if (highS < M and not is_gap(m_seq[x]))
        {
          highS = M;
          highX = x;
          highY = y;
        }
      }
      else if (Ix1 >= Iy1)
      {
        tb(x, y) = 1;
        B(x, y) = Ix1;

        Ix(x, y) = Ix1 - gep_a[x];
        Iy(x, y) = std::max(M - (y < dimY - 1 ? gop_b[y] : 0), Iy1 - gep_b[y]);
      }
      else
      {
        tb(x, y) = -1;
        B(x, y) = Iy1;

        Ix(x, y) = std::max(M - (x < dimX - 1 ? gop_a[x] : 0), Ix1 - gep_a[x]);
        Iy(x, y) = Iy1 - gep_b[y];
      }
    }
  }

//#if not NDEBUG
//
//bool dmp = false;
//if (e->m_acc == "Q0IJ18")
//{
//  dmp = true;
//  dump(B, Ix, Iy, tb, gop_a, gop_b, gep_a, gep_b, m_seq, e->m_seq);
//}
//
//#endif

  // build the alignment
  x = highX;
  y = highY;

  uint32 ident = 0;
  uint32 similar = 0;
  uint32 length = 0;
  uint32 lengthI = 0;
  uint32 xgaps = 0;
  uint32 xgapsI = 0;

  // trace back the matrix
  while (x >= 0 and y >= 0 and B(x, y) > 0)
  {
    ++lengthI;
    switch (tb(x, y))
    {
      case -1:
        --y;
        ++xgapsI;
        break;

      case 1:
        if (is_gap(m_seq[x]))
          --lengthI;
        --x;
        break;

      case 0:
        if (not is_gap(m_seq[x]))
        {
          length = lengthI;
          xgaps = xgapsI;

          if (e->m_seq[y] == m_seq[x])
            ++ident, ++similar;
          else if (score(kMPam250, m_seq[x], e->m_seq[y]) > 0)
            ++similar;
        }

        --x;
        --y;
        break;

      default:
        assert(false);
        break;
    }
  }

  uint32 tix = std::max(10U, std::min(length, 80U)) - 10;

  // Add the hit only if it is within the required parameters.
  //
  // accept only alignment long enough (suppress fragments)
  // and those that score high enough
  if (length >= m_seq.length() * m_frag_cutoff and
    ident >= length * (kHomologyThreshold[tix] + m_threshold))
  {
    // reserve space, if needed
    if (xgaps > 0)
    {
      const uint32 kBlockSize = 1024;
      uint32 n = static_cast<uint32>(
          (((m_residues.size() + xgaps) / kBlockSize) + 1) * kBlockSize);
      m_seq.reserve(n);
      foreach (MHitPtr e, m_entries)
        e->m_aligned.reserve(n);
    }

    int32 fx = x + 1, fy = y + 1;

    // update insert/delete counters for the residues
    x = highX;  e->m_ilas = m_residues[x].m_seq_nr;
    y = highY;  e->m_jlas = y + 1;

    // trace back to fill aligned sequence and to create gaps in MSA
    e->m_aligned = std::string(m_seq.length() + xgaps, '.');
    bool gappedx = false, gappedy = false;

    lengthI = length;
    while (x >= fx and y >= fy and lengthI-- > 0)
    {
      switch (tb(x, y))
      {
        case -1:
          e->m_aligned[x + xgaps] = kResidues[e->m_seq[y]];

          m_residues.insert(m_residues.begin() + x + 1,
            MResInfo::NewGap(m_entries.size() + 1,
                             m_sum_dist_weight, e->m_seq[y], e->m_distance));

          m_seq.insert(m_seq.begin() + x + 1, '.');

          foreach (MHitPtr e, m_entries)
            e->m_aligned.insert(e->m_aligned.begin() + x + 1, '.');

          --y;
          --xgaps;

          if (not gappedx)
            ++e->m_gaps;
          ++e->m_gapn;
          gappedx = true;
          gappedy = false;
          break;

        case 1:
          if (is_gap(m_seq[x]))
            ++lengthI;
          else
          {
            if (not gappedy)
              ++e->m_gaps;
            ++e->m_gapn;
            gappedx = false;
            gappedy = true;
          }

          m_residues[x].AddGap(e->m_distance);
          --x;
          break;

        case 0:
          e->m_aligned[x + xgaps] = kResidues[e->m_seq[y]];
          m_residues[x].Add(e->m_seq[y], e->m_distance);

          if (not is_gap(e->m_seq[y]) and is_gap(m_seq[x]))
          {
            if (not gappedx)
              ++e->m_gaps;
            gappedx = true;
            ++e->m_gapn;
          }
          else if (gappedx)
          {
            ++m_residues[x].m_ins;
            gappedx = false;
          }

          if (gappedy)
          {
            ++m_residues[x].m_del;
            gappedy = false;
          }

          --x;
          --y;
          break;
      }
    }

    // update the new entry
    e->m_identical = ident;
    e->m_similar = similar;
    e->m_length = length;
    e->m_distance = 1 - float(ident) / length;
    e->m_score = 1 - e->m_distance;

    e->m_ifir = m_residues[x + 1].m_seq_nr;
    e->m_jfir = y + 2;

    e->m_stid = e->m_acc + '/' +
      boost::lexical_cast<std::string>(e->m_jfir) + '-' +
      boost::lexical_cast<std::string>(e->m_jlas);

    m_entries.push_back(e);

    m_sum_dist_weight += e->m_distance;

//#if not defined(NDEBUG)
//    if (dmp)
//    {
//      std::string s = decode(m_seq);
//      for (std::string::size_type i = 72; i < s.length(); i += 73)
//        s.insert(s.begin() + i, '\n');
//
//      cout << '>' << "PDB" << std::endl
//         << s << std::endl;
//
//      foreach (MHitPtr e, m_entries)
//        cout << *e;
//      exit(1);
//    }
//#endif
  }
}

char map_value_to_char(uint32 v)
{
  char result = '0';
  if (v < 10)
    result += v;
  else
    result = '+';
  return result;
}

char map_value_to_char(double v)
{
  return map_value_to_char(static_cast<uint32>(v));
}

std::string FixedLengthString(const std::string &s, const uint64 length)
{
    std::string r(s);

    if (s.length() < length)
    {
        r = r.insert(0, length - r.length(), ' ');
    }
    else if (s.length() > length)
    {
        r = std::string(length, '-');
        r[length - 1] = '>';
    }
    return r;
}

std::string FixedLengthString(const int64 number, const uint64 length)
{
    std::string s = std::to_string(number);

    return FixedLengthString(s, length);
}

void MProfile::PrintStockholm(std::ostream& os, const std::string& inChainID,
                              bool inFetchDBRefs) const
{
  os << "#=GF ID " << inChainID << std::endl
     << "#=GF SQ " << m_entries.size() << std::endl;

  if (m_shuffled)
    os << "#=GF CC Since the number of hits exceeded the max-hits parameter,"
       << " a random set was chosen" << std::endl;

  // ## per residue information

  uint32 nextNr = m_residues.front().m_seq_nr;
  os << "#=GF CC ## RESIDUE INFORMATION" << std::endl
     << "#=GF CC SeqNo   PDBNo AA STRUCTURE BP1 BP2  ACC  NOCC VAR CHAIN AUTHCHAIN     NUMBER     RESNUM        BP1        BP2"
     << std::endl;
  for (const MResInfo &ri : m_residues)
  {
    if (ri.m_chain_id.empty())
      continue;

    if (ri.m_seq_nr != nextNr)
      os << boost::format("#=GF RI %5.5d       ! !              0   0    0     0   0") % nextNr << std::endl;

    uint32 ivar = uint32(100 * (1 - ri.m_consweight));

    char aa_code = GetAACode(ri.m_aa, ri.m_ss_bridge_nr),
         ins_code = (ri.m_ins_code.size() > 0)? ri.m_ins_code.at(0): ' ';

    os << boost::format("#=GF RI %5.5s %5.5s%c%c %c %9.9s %4.4s%4.4s%5.5s %5.5d%4.4d")
          % FixedLengthString(ri.m_seq_nr, 5) % FixedLengthString(ri.m_pdb_nr, 5) % ins_code % FixedLengthString(ri.m_chain_id, 1)
          % aa_code % ri.m_structure
          % FixedLengthString(ri.m_beta_partner_1, 4) % FixedLengthString(ri.m_beta_partner_2, 4)
          % FixedLengthString(floor(ri.m_accessibility + 0.5), 5)
          % ri.m_nocc % ivar;

    os << boost::format("  %4.4s      %4.4s %10.10s %10.10s %10.10s %10.10s")
          % ri.m_chain_id % ri.m_auth_chain_id
          % FixedLengthString(ri.m_seq_nr, 10) % FixedLengthString(ri.m_pdb_nr, 10)
          % FixedLengthString(ri.m_beta_partner_1, 10) % FixedLengthString(ri.m_beta_partner_2, 10) << std::endl;

    nextNr = ri.m_seq_nr + 1;
  }

  // ## SEQUENCE PROFILE AND ENTROPY
  os << "#=GF CC ## SEQUENCE PROFILE AND ENTROPY" << std::endl
     << "#=GF CC   SeqNo PDBNo   V   L   I   M   F   W   Y   G   A   P   S   T   C   H   R   K   Q   E   N   D  NOCC NDEL NINS ENTROPY RELENT WEIGHT CHAIN AUTHCHAIN     NUMBER     RESNUM        BP1        BP2" << std::endl;

  nextNr = m_residues.front().m_seq_nr;
  for (const MResInfo &ri : m_residues)
  {
    if (ri.m_chain_id.empty())
      continue;

    if (ri.m_seq_nr != nextNr)
      os << boost::format("#=GF PR %5.5d           0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0     0    0    0   0.000      0  1.00")
        % nextNr << std::endl;

    char chainChar = '>';
    if (ri.m_chain_id.length() <= 1)
        chainChar = ri.m_chain_id[0];

    os << boost::format("#=GF PR %5.5d %5.5d%1.1s%1.1s") % ri.m_seq_nr % ri.m_pdb_nr % ri.m_ins_code % chainChar;

    for (uint32 i = 0; i < 20; ++i)
      os << boost::format("%4.4d") % uint32(100.0 * ri.m_freq[i] + 0.5);

    uint32 relent = uint32(100 * ri.m_entropy / log(20.0));
    os << "  " << boost::format("%4.4d %4.4d %4.4d  %6.3f   %4.4d %5.2f") % ri.m_nocc % ri.m_del % ri.m_ins % ri.m_entropy % relent % ri.m_consweight;

    os << boost::format("   %4.4s      %4.4s %10.10s %10.10s %10.10s %10.10s")
          % ri.m_chain_id % ri.m_auth_chain_id
          % FixedLengthString(ri.m_seq_nr, 10) % FixedLengthString(ri.m_pdb_nr, 10)
          % FixedLengthString(ri.m_beta_partner_1, 10) % FixedLengthString(ri.m_beta_partner_2, 10) << std::endl;

    nextNr = ri.m_seq_nr + 1;
  }

  // find the longest ID std::string length
  std::string::size_type tl = inChainID.length();
  foreach (const MHitPtr e, m_entries)
  {
    if (tl < e->m_stid.length())
      tl = e->m_stid.length();
  }

  os << "#=GS " << inChainID << std::string(tl - inChainID.length(), ' ') << " CC The query chain" << std::endl;

  boost::format fmt("#=GS %s HSSP score=%4.2f/%4.2f aligned=%d-%d/%d-%d length=%d ngaps=%d gaplen=%d seqlen=%d");

  std::map<std::string,std::vector<std::string>> linked;
#ifdef HAVE_LIBZEEP
  if (inFetchDBRefs)
  {
    const std::string kBaseURL = "http://mrs.cmbi.umcn.nl/mrsws/search";

    foreach (const MHitPtr e, m_entries)
      linked[e->m_id].clear();

    FetchPDBReferences(kBaseURL, "uniprot", linked);
  }
#endif

  foreach (const MHitPtr e, m_entries)
  {
    std::string id = e->m_stid + std::string(tl - e->m_stid.length(), ' ');

    os << "#=GS " << id << " ID " << e->m_id << std::endl
       << "#=GS " << id << " DE " << e->m_def << std::endl
       << fmt % id % e->m_score % (float(e->m_similar) / e->m_length)
           % e->m_ifir % e->m_ilas % e->m_jfir % e->m_jlas % e->m_length
           % e->m_gaps % e->m_gapn % e->m_seq.length() << std::endl;

#ifdef HAVE_LIBZEEP
    if (inFetchDBRefs and not linked[e->m_id].empty())
    {
      os << "#=GS " << id << " DR PDB " << ba::join(linked[e->m_id], ", ")
         << std::endl;
    }
#endif
  }

  if (tl < 17)
    tl = 17;

  std::string::size_type o = 0;
  while (o < m_seq.length())
  {
    std::string::size_type n = 72;
    if (o + n > m_seq.length())
      n = m_seq.length() - o;

    os << std::endl
       << inChainID << std::string(tl - inChainID.length() + 1, ' ')
       << decode(m_seq.substr(o, n)) << std::endl;

    std::string ss(n, '.'), ins(n, ' '), del(n, ' '), ent(n, '-'), var(n, '-');
    for (std::string::size_type i = o; i < o + n; ++i)
    {
      if (m_residues[i].m_seq_nr == 0)
        continue;

      switch (m_residues[i].m_ss)
      {
        case alphahelix: ss[i - o] = 'H'; break;
        case betabridge: ss[i - o] = 'B'; break;
        case strand:     ss[i - o] = 'E'; break;
        case helix_3:    ss[i - o] = 'G'; break;
        case helix_5:    ss[i - o] = 'I'; break;
        case turn:       ss[i - o] = 'T'; break;
        case bend:       ss[i - o] = 'S'; break;
        case loop:       ss[i - o] = 'C'; break;
      }

      ent[i - o] = map_value_to_char(10 * m_residues[i].m_entropy / log(20.0));
      var[i - o] = map_value_to_char(10 * (1 - m_residues[i].m_consweight));
    }

    foreach (const MHitPtr e, m_entries)
      os << e->m_stid << std::string(tl - e->m_stid.length() + 1, ' ')
         << e->m_aligned.substr(o, n) << std::endl;

    os << "#=GC SS          " << std::string(tl - 17 + 1, ' ') << ss
       << std::endl
       << "#=GC Entropy     " << std::string(tl - 17 + 1, ' ') << ent
       << std::endl
       << "#=GC Variability " << std::string(tl - 17 + 1, ' ') << var
       << std::endl;

    o += n;
  }

  os << "//" << std::endl;
}

void MProfile::PrintStockholm(std::ostream& os, const MProtein& inProtein,
                              bool inFetchDBRefs,
                              const std::vector<std::string>& inUsed,
                              const std::vector<std::string>& inAKA) const
{
  using namespace boost::gregorian;
  date today = day_clock::local_day();

  // write out the profile in Stockholm 1.0 format

  os << "# STOCKHOLM 1.0" << std::endl
     << "#=GF CC DATE   " << to_iso_extended_string(today) << std::endl;

  std::string s = inProtein.GetID();
  if (not s.empty())
    os << "#=GF CC PDBID  " << s << std::endl;

  s = inProtein.GetHeader();
  if (not s.empty())
    os << "#=GF CC HEADER " << s.substr(10) << std::endl;

  s = inProtein.GetCompound();
  if (not s.empty())
    os << "#=GF CC COMPND " << s.substr(10) << std::endl;

  s = inProtein.GetSource();
  if (not s.empty())
    os << "#=GF CC SOURCE " << s.substr(10) << std::endl;

  s = inProtein.GetAuthor();
  if (not s.empty())
    os << "#=GF CC AUTHOR " << s.substr(10) << std::endl;

  foreach (auto dbref, inProtein.GetDbRef())
    os << "#=GF CC " << dbref << std::endl;

  std::string queryID = inProtein.GetID();
  if (inProtein.GetChains().size() > 1)
  {
    auto fmt = [](const std::vector<std::string>& a) -> std::string {
      std::string result;
      if (a.size() == 1)
        result += a.front();
      else if (not a.empty())
      {
        result = ba::join(a | boost::adaptors::sliced(0, a.size() - 1), ", ");
        result += " and ";
        result += a.back();
      }
      return result;
    };

    os << "#=GF CC PDB file contains " << inProtein.GetChains().size()
       << " chains. Used chain" << ( inUsed.size() > 1 ? "s are " : " is " )
       << fmt(inUsed) << std::endl;
    if (not inAKA.empty())
      os << "#=GF CC Chain " << m_chain.GetChainID()
         << " is considered to be the same as " << fmt(inAKA) << std::endl;
    queryID = inProtein.GetID() + '/' + m_chain.GetChainID();
  }

  PrintStockholm(os, queryID, inFetchDBRefs);
}

// --------------------------------------------------------------------

void MProfile::Process(std::istream& inHits, float inGapOpen,
                       float inGapExtend, uint32 inMaxHits, uint32 inThreads)
{
  std::vector<MHitPtr> hits;

  std::string id, def, seq;
  for (;;)
  {
    std::string line;
    getline(inHits, line);
    if (line.empty() and inHits.eof())
      break;

    if (ba::starts_with(line, ">"))
    {
      if (not (id.empty() or seq.empty()))
        hits.push_back(MHit::Create(id, def, seq));

      id.clear();
      def.clear();
      seq.clear();

      std::string::size_type s = line.find(' ');
      if (s != std::string::npos)
      {
        id = line.substr(1, s - 1);
        def = line.substr(s + 1);
      }
      else
        id = line.substr(1);
    }
    else
      seq += line;
  }

  if (not (id.empty() or seq.empty()))
    hits.push_back(MHit::Create(id, def, seq));

  // Now calculate distances
  MProgress p1(hits.size(), "distance");

  boost::thread_group threads;
  MCounter ix(0);

  for (uint32 t = 0; t < inThreads; ++t)
    threads.create_thread([this, &ix, &hits, &p1]() {
      for (;;)
      {
        uint64 next = ix++;
        if (next >= hits.size())
          break;

        hits[next]->CalculateDistance(m_seq);
        p1.Consumed(1);
      }
    });

  threads.join_all();

  // if we have way too many hits, take a random set
  if (hits.size() > inMaxHits * 10 and inMaxHits > 0)
  {
    if (VERBOSE)
      std::cerr << "dropping " << (hits.size() - 10 * inMaxHits) << " hits"
                << std::endl;

    random_shuffle(hits.begin(), hits.end());
    hits.erase(hits.begin() + inMaxHits * 10, hits.end());
    m_shuffled = true;
  }

  // sort them by distance
  sort(hits.begin(), hits.end(), [](const MHitPtr a, const MHitPtr b) -> bool {
    return a->m_distance < b->m_distance;
  });

  // and then align all the hits
  MProgress p2(hits.size(), "aligning");
  foreach (MHitPtr e, hits)
  {
    Align(e, inGapOpen, inGapExtend);
    p2.Consumed(1);
  }

  // now if we have too many entries, take a random set
  if (m_entries.size() > inMaxHits and inMaxHits > 0)
  {
    random_shuffle(m_entries.begin(), m_entries.end());
    m_entries.erase(m_entries.begin() + inMaxHits, m_entries.end());
    m_shuffled = true;
  }

  // sort by score
  sort(m_entries.begin(), m_entries.end(),
       [](const MHitPtr a, const MHitPtr b) -> bool {
    return a->m_score > b->m_score;
  });

  CalculateConservation(inThreads);
}

// --------------------------------------------------------------------

// Find the minimal set of overlapping sequences
// In case of strong similarity (distance <= 0.01) we take the longest chain.
void ClusterSequences(const std::vector<sequence>& s, std::vector<size_t>& ix)
{
  std::vector<bool> skip(s.size(), false);

  for (;;)
  {
    bool found = false;
    for (uint32 i = 0; not found and i < s.size() - 1; ++i)
    {
      for (uint32 j = i + 1; not found and j < s.size(); ++j)
      {
        if (skip[i] or skip[j])
          continue;

        const sequence& a = s[i];
        const sequence& b = s[j];

        bool isSame = a == b;
        if (not isSame)
        {
          float d = calculateDistance(a, b);
          // rescale distance to shortest length:
          d = 1 - (1 - d) * std::max(a.length(), b.length()) / std::min(a.length(), b.length());
          isSame = (d <= 0.01);
        }

        if (isSame)
        {
          skip[j] = true;
          ix[j] = i;
          found = true;
        }
      }
    }

    if (not found)
      break;
  }

  // change ix entries to make sure we use only the longest chains
  for (uint32 i = 0; i < ix.size(); ++i)
  {
    if (ix[i] == i)
      continue;

    uint32 m = ix[i];

    uint32 la = s[m].length();
    uint32 lb = s[i].length();
    if (la < lb)
    {
      for_each(ix.begin(), ix.end(), [=](size_t& j) { if (j == m) j = i; });
    }
  }
}

// --------------------------------------------------------------------
// Calculate the variability of a residue, based on dayhoff similarity
// and weights

std::pair<const char*,uint32> kSentinel((const char*)nullptr, 0);

void CalculateConservation(buffer<std::pair<const char*,uint32>>& b,
                           const std::vector<MHitPtr>& inHits,
                           std::vector<float>& sumvar,
                           std::vector<float>& sumdist)
{
  size_t length = sumvar.size();
  std::vector<float> simval(length);

  for (;;)
  {
    auto next = b.get();
    if (next == kSentinel)
      break;

    const char* si = next.first;
    for (uint32 j = next.second; j < inHits.size(); ++j)
    {
      const char* sj = inHits[j]->m_aligned.c_str();

      uint32 len = 0, agr = 0;
      for (uint32 k = 0; k < length; ++k)
      {
        simval[k] = std::numeric_limits<float>::min();

        if (is_gap(si[k]) or is_gap(sj[k]))
          continue;

        ++len;
        if (si[k] == sj[k])
          ++agr;

        int8 ri = ResidueNr(si[k]);
        int8 rj = ResidueNr(sj[k]);

        if (ri <= 20 and rj <= 20)
          simval[k] = score(kDayhoffData, ri, rj);
      }

      if (len == 0)
        continue;

      float distance = 1 - (float(agr) / float(len));
      for (uint32 k = 0; k < length; ++k)
      {
        if (simval[k] != std::numeric_limits<float>::min())
        {
          sumvar[k] += distance * simval[k];
          sumdist[k] += distance * 1.5f;
        }
      }
    }
  }

  b.put(kSentinel);
}

void MProfile::CalculateConservation(uint32 inThreads)
{
  std::vector<float> sumvar(m_seq.length(), 0), sumdist(m_seq.length(), 0);

  // Calculate conservation weights in multiple threads to gain speed.
  buffer<std::pair<const char*,uint32>> b;
  boost::thread_group threads;
  boost::mutex sumLock;

  for (uint32 t = 0; t < inThreads; ++t)
    threads.create_thread([&]() {
      std::vector<float> csumvar(sumvar.size(), 0), csumdist(sumdist.size(), 0);

      HSSP::CalculateConservation(b, m_entries, csumvar, csumdist);

      // accumulate our data
      boost::mutex::scoped_lock l(sumLock);

      for (size_t i = 0; i < sumvar.size(); ++i)
      {
        sumvar[i] += csumvar[i];
        sumdist[i] += csumdist[i];
      }
    });

  int64 N = (m_entries.size() * (m_entries.size() + 1)) / 2;
  if (m_shuffled)
    N += m_seq.length();

  MProgress p(N, "conservation");

  if (m_shuffled)  // need to recalculate m_dist[] and m_nocc
  {
    for (uint32 i = 0; i < m_seq.length(); ++i)
    {
      MResInfo& ri = m_residues[i];

      ri.m_nocc = is_gap(m_seq[i]) ? 0 : 1;
      for (uint32 j = 0; j < 23; ++j)
        ri.m_dist[j] = m_seq[i] == j ? 1 : 0;

      foreach (MHitPtr e, m_entries)
      {
        uint8 r = ResidueNr(e->m_aligned[i]);
        if (r < 23)
        {
          ++ri.m_dist[r];
          ++ri.m_nocc;
        }
      }

      p.Consumed(1);
    }
  }

  std::string s(decode(m_seq));
  b.put(std::make_pair(s.c_str(), 0));

  p.Consumed(m_entries.size());

  for (uint32 i = 0; i + 1 < m_entries.size(); ++i)
  {
    b.put(std::make_pair(m_entries[i]->m_aligned.c_str(), i + 1));
    p.Consumed(m_entries.size() - i);
  }

  b.put(kSentinel);
  threads.join_all();

  for (uint32 i = 0; i < m_seq.length(); ++i)
  {
    MResInfo& ri = m_residues[i];

    if (sumdist[i] > 0)
      ri.m_consweight = std::min(1.0f, sumvar[i] / sumdist[i]);
    else
      ri.m_consweight = 1;

    ri.m_entropy = 0;

    static const int8 kResIx[] = {
      //  V   L   I   M   F   W   Y   G   A   P   S   T   C   H   R   K   Q   E   N   D
//         18, 10,  8, 11,  5, 19, 20,  6,  0, 13, 16, 17,  2,  7, 15,  9, 14,  4, 12,  3
         17,  9,  7, 10,  4, 18, 19,  5,  0, 12, 15, 16,  1,  6, 14,  8, 13,  3, 11,  2
    };

    for (uint32 i = 0; i < 20; ++i)
    {
      ri.m_freq[i] = float(ri.m_dist[kResIx[i]]) / ri.m_nocc;
      if (ri.m_freq[i] > 0)
        ri.m_entropy -= ri.m_freq[i] * log(ri.m_freq[i]);
    }
  }
}

// --------------------------------------------------------------------

void CreateHSSP(const MProtein& inProtein,
                const std::vector<fs::path>& inDatabanks,
                uint32 inMaxHits, uint32 inMinSeqLength, float inGapOpen,
                float inGapExtend, float inThreshold, float inFragmentCutOff,
                uint32 inThreads, bool inFetchDBRefs, std::ostream& inOs)
{
  // construct a set of unique sequences, containing only the largest ones in
  // case of overlap
  std::vector<sequence> seqset;
  std::vector<size_t> ix;
  std::vector<const MChain*> chains;
  std::vector<std::vector<std::string>> aka;
  std::vector<std::string> used;

  foreach (const MChain* chain, inProtein.GetChains())
  {
    std::string seq;
    chain->GetSequence(seq);

    if (seq.length() < inMinSeqLength)
      continue;

    chains.push_back(chain);
    seqset.push_back(encode(seq));
    ix.push_back(ix.size());
    aka.push_back(std::vector<std::string>());
  }

  if (seqset.empty())
    throw mas_exception(boost::format("Not enough sequences in PDB file of length %1%") % inMinSeqLength);

  if (seqset.size() > 1)
    ClusterSequences(seqset, ix);

  // only take the unique sequences
  for (size_t i = 0; i < ix.size(); ++i)
  {
    if (ix[i] != i)
      aka[ix[i]].push_back(chains[i]->GetChainID());
  }

  sort(ix.begin(), ix.end());
  ix.erase(unique(ix.begin(), ix.end()), ix.end());

  foreach (size_t i, ix)
    used.push_back(chains[i]->GetChainID());

  bool empty = true;

  foreach (size_t i, ix)
  {
    const MChain& chain(*chains[i]);

    // do a blast search for inMaxHits * 4 hits.
    std::vector<char> blastHits;

    io::filtering_ostream out(io::back_inserter(blastHits));

    //ifstream f("1F88-A-hits.fa");
    //io::copy(f, out);

    std::string seq = decode(seqset[i]);
    SearchAndWriteResultsAsFastA(out, inDatabanks, seq,
      "blastp", "BLOSUM62", 3, 10, true, true, -1, -1, 0, inThreads);
    out.flush();

    if (blastHits.empty())
      continue;

    MProfile profile(chain, seqset[i], inThreshold, inFragmentCutOff);

    io::filtering_istream in(boost::make_iterator_range(blastHits));
    profile.Process(in, inGapOpen, inGapExtend, inMaxHits, inThreads);

    if (profile.m_entries.empty())
      continue;

    empty = false;
    profile.PrintStockholm(inOs, inProtein, inFetchDBRefs, used, aka[i]);
  }

  if (empty)
    throw mas_exception("No hits found");
}

// --------------------------------------------------------------------

void CreateHSSP(const std::string& inProtein,
                const std::vector<fs::path>& inDatabanks,
                uint32 inMaxHits,
                uint32 inMinSeqLength, float inGapOpen, float inGapExtend,
                float inThreshold, float inFragmentCutOff, uint32 inThreads,
                bool inFetchDBRefs, std::ostream& inOs)
{
  MChain* chain = new MChain("A");
  std::vector<MResidue*>& residues = chain->GetResidues();
  MResidue* last = nullptr;
  int32 nr = 1;
  foreach (char r, inProtein)
  {
    residues.push_back(new MResidue(nr, r, last));
    ++nr;
    last = residues.back();
  }

  MProtein protein("INPUT", chain);
  CreateHSSP(protein, inDatabanks, inMaxHits, inMinSeqLength, inGapOpen,
             inGapExtend, inThreshold, inFragmentCutOff, inThreads,
             inFetchDBRefs, inOs);
}

}

