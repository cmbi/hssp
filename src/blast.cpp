// Copyright Maarten L. Hekkelman, Radboud University 2008-2011.
// Copyright Coos Baakman, Jon Black, Wouter G. Touw & Gert Vriend, Radboud university medical center 2015.
//   Distributed under the Boost Software License, Version 1.0.
//       (See accompanying file LICENSE_1_0.txt or copy at
//             http://www.boost.org/LICENSE_1_0.txt)

#include "blast.h"

#include "matrix.h"
#include "utils.h"
#include "progress.h"

#include <boost/algorithm/string.hpp>
#include <boost/detail/atomic_count.hpp>
#include <boost/foreach.hpp>
#include <boost/format.hpp>
#include <boost/iostreams/device/mapped_file.hpp>
#include <boost/regex.hpp>

#include <limits>
#include <numeric>

namespace ba = boost::algorithm;
namespace fs = boost::filesystem;
namespace io = boost::iostreams;


#define foreach BOOST_FOREACH
// --------------------------------------------------------------------

boost::regex
  kFastARE("^>(\\w+)((?:\\|([^| ]*))?(?:\\|([^| ]+))?(?:\\|([^| ]+))?(?:\\|([^| ]+))?)(?: (.+))\n?");

const uint32
  kAACount        = 22,  // 20 + B and Z
  kResCount        = 23,  // includes X
  kBits          = 5,
  kThreshold        = 11,
  kUngappedDropOff    = 7,
  kGappedDropOff      = 15,
  kGappedDropOffFinal    = 25,
  kGapTrigger        = 22,

  kMaxSequenceLength    = std::numeric_limits<uint16>::max();

const int32
  kHitWindow        = 40;

const double
  kLn2          = log(2.);

const int16
  kSentinalScore      = -9999;

class Matrix
{
  public:
    Matrix(const std::string& inName, int32 inGapOpen, int32 inGapExtend);

  int8 operator()(char inAA1, char inAA2) const;
  int8 operator()(uint8 inAA1, uint8 inAA2) const;

  int32 OpenCost() const { return mData.mGapOpen; }
  int32 ExtendCost() const { return mData.mGapExtend; }

  double GappedLambda() const { return mData.mGappedStats.lambda; }
  double GappedKappa() const { return mData.mGappedStats.kappa; }
  double GappedEntropy() const { return mData.mGappedStats.entropy; }
  double GappedAlpha() const { return mData.mGappedStats.alpha; }
  double GappedBeta() const { return mData.mGappedStats.beta; }

  double UngappedLambda() const { return mData.mUngappedStats.lambda; }
  double UngappedKappa() const { return mData.mUngappedStats.kappa; }
  double UngappedEntropy() const { return mData.mUngappedStats.entropy; }
  double UngappedAlpha() const { return mData.mUngappedStats.alpha; }
  double UngappedBeta() const { return mData.mUngappedStats.beta; }

  private:
    MMatrixData mData;
};

Matrix::Matrix(const std::string& inName, int32 inGapOpen, int32 inGapExtend)
{
  mData.mName = nullptr;
  for (const MMatrixData* data = kMMatrixData; data->mName != nullptr; ++data)
  {
    if (ba::iequals(inName, data->mName) and
      inGapOpen == data->mGapOpen and
      inGapExtend == data->mGapExtend)
    {
      mData = *data;
      break;
    }
  }

  if (mData.mName == nullptr)
    throw mas_exception(boost::format("Unsupported matrix/gap combination (%s/%d/%d)") % inName % inGapOpen % inGapExtend);
}

inline int8 Matrix::operator()(uint8 inAA1, uint8 inAA2) const
{
  int8 result;

  if (inAA1 >= inAA2)
    result = mData.mMatrix[(inAA1 * (inAA1 + 1)) / 2 + inAA2];
  else
    result = mData.mMatrix[(inAA2 * (inAA2 + 1)) / 2 + inAA1];

  return result;
}

inline int8 Matrix::operator()(char inAA1, char inAA2) const
{
  return operator()(ResidueNr(inAA1), ResidueNr(inAA2));
}

// --------------------------------------------------------------------

namespace filter
{

class Alphabet
{
  public:
    Alphabet(const char* inChars);

  bool Contains(char inChar) const;
  long GetIndex(char inChar) const;
  long GetSize() const { return mAlphaSize; }
  double GetLnSize() const { return mAlphaLnSize; }

  private:
  long mAlphaSize;
  double mAlphaLnSize;
  long mAlphaIndex[128];
  const char* mAlphaChars;
};

Alphabet::Alphabet(const char* inChars)
  : mAlphaChars(inChars)
{
  mAlphaSize = static_cast<int>(strlen(inChars));
  mAlphaLnSize = log(static_cast<double>(mAlphaSize));

  for (uint32 i = 0; i < 128; ++i)
  {
    mAlphaIndex[i] = static_cast<long>(
        std::find(mAlphaChars, mAlphaChars + mAlphaSize,
                  toupper(i)) - mAlphaChars);
  }
}

bool Alphabet::Contains(char inChar) const
{
  bool result = false;
  if (inChar >= 0)
    result = mAlphaIndex[toupper(inChar)] < mAlphaSize;
  return result;
}

long Alphabet::GetIndex(char inChar) const
{
  return mAlphaIndex[toupper(inChar)];
}

const Alphabet
  kProtAlphabet = Alphabet("ACDEFGHIKLMNPQRSTVWY"),
  kNuclAlphabet = Alphabet("ACGTU");

class Window
{
  public:
      Window(const std::string& inSequence, long inStart, long inLength,
             const Alphabet& inAlphabet);

  void CalcEntropy();
  bool ShiftWindow();

  double GetEntropy() const { return mEntropy; }
  long GetBogus() const { return mBogus; }

  void DecState(long inCount);
  void IncState(long inCount);

  void Trim(long& ioEndL, long& ioEndR, long inMaxTrim);

  private:
  const std::string& mSequence;
  std::vector<long> mComposition;
  std::vector<long> mState;
  long mStart;
  long mLength;
  long mBogus;
  double mEntropy;
  const Alphabet&mAlphabet;
};

Window::Window(const std::string& inSequence, long inStart, long inLength,
               const Alphabet& inAlphabet)
  : mSequence(inSequence),
    mComposition(inAlphabet.GetSize()),
    mStart(inStart),
    mLength(inLength),
    mBogus(0),
    mEntropy(-2.0),
    mAlphabet(inAlphabet)
{
  long alphaSize = mAlphabet.GetSize();

  for (long i = mStart; i < mStart + mLength; ++i)
  {
    if (mAlphabet.Contains(mSequence[i]))
      ++mComposition[mAlphabet.GetIndex(mSequence[i])];
    else
      ++mBogus;
  }

  mState.insert(mState.begin(), alphaSize + 1, 0);

  int n = 0;
  for (long i = 0; i < alphaSize; ++i)
  {
    if (mComposition[i] > 0)
    {
      mState[n] = mComposition[i];
      ++n;
    }
  }

  std::sort(mState.begin(), mState.begin() + n, std::greater<long>());
}

void Window::CalcEntropy()
{
  mEntropy = 0.0;

  double total = 0.0;
  for (uint32 i = 0; i < mState.size() and mState[i] != 0; ++i)
    total += mState[i];

  if (total != 0.0)
  {
    for (uint32 i = 0; i < mState.size() and mState[i]; ++i)
    {
      double t = mState[i] / total;
      mEntropy += t * log(t);
    }
    mEntropy = fabs(mEntropy / log(0.5));
  }
}

void Window::DecState(long inClass)
{
  for (uint32 ix = 0; ix < mState.size() and mState[ix] != 0; ++ix)
  {
    if (mState[ix] == inClass and mState[ix + 1] < inClass)
    {
      --mState[ix];
      break;
    }
  }
}

void Window::IncState(long inClass)
{
  for (uint32 ix = 0; ix < mState.size(); ++ix)
  {
    if (mState[ix] == inClass)
    {
      ++mState[ix];
      break;
    }
  }
}

bool Window::ShiftWindow()
{
  if (uint32(mStart + mLength) >= mSequence.length())
    return false;

  char ch = mSequence[mStart];
  if (mAlphabet.Contains(ch))
  {
    long ix = mAlphabet.GetIndex(ch);
    DecState(mComposition[ix]);
    --mComposition[ix];
  }
  else
    --mBogus;

  ++mStart;

  ch = mSequence[mStart + mLength - 1];
  if (mAlphabet.Contains(ch))
  {
    long ix = mAlphabet.GetIndex(ch);
    IncState(mComposition[ix]);
    ++mComposition[ix];
  }
  else
    ++mBogus;

  if (mEntropy > -2.0)
    CalcEntropy();

  return true;
}

static double lnfac(long inN)
{
    const double c[] = {
         76.18009172947146,
        -86.50532032941677,
         24.01409824083091,
        -1.231739572450155,
         0.1208650973866179e-2,
        -0.5395239384953e-5
    };
  static std::map<long,double> sLnFacMap;

  if (sLnFacMap.find(inN) == sLnFacMap.end())
  {
    double x = inN + 1;
    double t = x + 5.5;
    t -= (x + 0.5) * log(t);
      double ser = 1.000000000190015;
      for (int i = 0; i <= 5; i++)
      {
        ++x;
          ser += c[i] / x;
      }
      sLnFacMap[inN] = -t + log(2.5066282746310005 * ser / (inN + 1));
  }

  return sLnFacMap[inN];
}

static double lnperm(std::vector<long>& inState, long inTotal)
{
  double ans = lnfac(inTotal);
  for (uint32 i = 0; i < inState.size() and inState[i] != 0; ++i)
    ans -= lnfac(inState[i]);
  return ans;
}

static double lnass(std::vector<long>& inState, Alphabet inAlphabet)
{
    double result = lnfac(inAlphabet.GetSize());
    if (inState.size() == 0 or inState[0] == 0)
        return result;

    int total = inAlphabet.GetSize();
    int cl = 1;
    int i = 1;
    int sv_cl = inState[0];

    while (inState[i] != 0)
    {
        if (inState[i] == sv_cl)
            cl++;
        else
        {
            total -= cl;
            result -= lnfac(cl);
            sv_cl = inState[i];
            cl = 1;
        }
        i++;
    }

    result -= lnfac(cl);
    total -= cl;
    if (total > 0)
      result -= lnfac(total);

    return result;
}

static double lnprob(std::vector<long>& inState, long inTotal,
                     const Alphabet& inAlphabet)
{
  double ans1, ans2 = 0, totseq;

  totseq = inTotal * inAlphabet.GetLnSize();
  ans1 = lnass(inState, inAlphabet);
  if (ans1 > -100000.0 and inState[0] != std::numeric_limits<long>::min())
    ans2 = lnperm(inState, inTotal);
  else
    throw mas_exception("Error in calculating lnass");
  return ans1 + ans2 - totseq;
}

void Window::Trim(long& ioEndL, long& ioEndR, long inMaxTrim)
{
  double minprob = 1.0;
  long lEnd = 0;
  long rEnd = mLength - 1;
  int minLen = 1;
  int maxTrim = inMaxTrim;
  if (minLen < mLength - maxTrim)
    minLen = mLength - maxTrim;

  for (long len = mLength; len > minLen; --len)
  {
    Window w(mSequence, mStart, len, mAlphabet);

    int i = 0;
    bool shift = true;
    while (shift)
    {
      double prob = lnprob(w.mState, len, mAlphabet);
      if (prob < minprob)
      {
        minprob = prob;
        lEnd = i;
        rEnd = len + i - 1;
      }
      shift = w.ShiftWindow();
      ++i;
    }
  }

  ioEndL += lEnd;
  ioEndR -= mLength - rEnd - 1;
}

static bool GetEntropy(const std::string& inSequence,
                       const Alphabet& inAlphabet, long inWindow,
                       long inMaxBogus, std::vector<double>& outEntropy)
{
  bool result = false;

  long downset = (inWindow + 1) / 2 - 1;
  long upset = inWindow - downset;

  if (inWindow <= static_cast<long>(inSequence.length()))
  {
    result = true;
    outEntropy.clear();
    outEntropy.insert(outEntropy.begin(), inSequence.length(), -1.0);

    Window win(inSequence, 0, inWindow, inAlphabet);
    win.CalcEntropy();

    long first = downset;
    long last = static_cast<long>(inSequence.length() - upset);
    for (long i = first; i <= last; ++i)
    {
//      if (GetPunctuation() and win.HasDash())
//      {
//        win.ShiftWindow();
//        continue;
//      }
      if (win.GetBogus() > inMaxBogus)
        continue;

      outEntropy[i] = win.GetEntropy();
      win.ShiftWindow();
    }
  }

  return result;
}

static void GetMaskSegments(bool inProtein, const std::string& inSequence,
                            long inOffset,
                            std::vector<std::pair<long,long> >& outSegments)
{
  double loCut, hiCut;
  long window, maxbogus, maxtrim;
  const Alphabet* alphabet;

  if (inProtein)
  {
    window = 12;
    loCut = 2.2;
    hiCut = 2.5;
    maxtrim = 50;
    maxbogus = 2;
    alphabet = &kProtAlphabet;
  }
  else
  {
    window = 32;
    loCut = 1.4;
    hiCut = 1.6;
    maxtrim = 100;
    maxbogus = 3;
    alphabet = &kNuclAlphabet;
  }

  long downset = (window + 1) / 2 - 1;
  long upset = window - downset;

  std::vector<double> e;
  GetEntropy(inSequence, *alphabet, window, maxbogus, e);

  long first = downset;
  long last = static_cast<long>(inSequence.length() - upset);
  long lowlim = first;

  for (long i = first; i <= last; ++i)
  {
    if (e[i] <= loCut and e[i] != -1.0)
    {
      long loi = i;
      while (loi >= lowlim and e[loi] != -1.0 and e[loi] <= hiCut)
        --loi;
      ++loi;

      long hii = i;
      while (hii <= last and e[hii] != -1.0 and e[hii] <= hiCut)
        ++hii;
      --hii;

      long leftend = loi - downset;
      long rightend = hii + upset - 1;

      std::string s(inSequence.substr(leftend, rightend - leftend + 1));
      Window w(s, 0, rightend - leftend + 1, *alphabet);
      w.Trim(leftend, rightend, maxtrim);

      if (i + upset - 1 < leftend)
      {
        long lend = loi - downset;
        long rend = leftend - 1;

        std::string left(inSequence.substr(lend, rend - lend + 1));
        GetMaskSegments(inProtein, left, inOffset + lend, outSegments);
      }

      outSegments.push_back(
        std::pair<long,long>(leftend + inOffset, rightend + inOffset + 1));
      i = rightend + downset;
      if (i > hii)
        i = hii;
      lowlim = i + 1;
    }
  }
}

std::string SEG(const std::string& inSequence)
{
  std::string result = inSequence;

  std::vector<std::pair<long,long> > segments;
  GetMaskSegments(true, result, 0, segments);

  for (uint32 i = 0; i < segments.size(); ++i)
  {
    for (long j = segments[i].first; j < segments[i].second; ++j)
      result[j] = 'X';
  }

  return result;
}

std::string DUST(const std::string& inSequence)
{
  std::string result = inSequence;

  std::vector<std::pair<long,long> > segments;
  GetMaskSegments(false, inSequence, 0, segments);

  for (uint32 i = 0; i < segments.size(); ++i)
  {
    for (long j = segments[i].first; j < segments[i].second; ++j)
      result[j] = 'X';
  }

  return result;
}

//int main()
//{
//  std::string seq;
//
//  ifstream in("input.seq", ios::binary);
//  in >> seq;
//  cout << FilterProtSeq(seq);
//  return 0;
//}


}

// --------------------------------------------------------------------

namespace ncbi
{

/**
 * Computes the adjustment to the lengths of the query and database sequences
 * that is used to compensate for edge effects when computing evalues.
 *
 * The length adjustment is an integer-valued approximation to the fixed
 * point of the function
 *
 *    f(ell) = beta +
 *               (alpha/lambda) * (log K + log((m - ell)*(n - N ell)))
 *
 * where m is the query length n is the length of the database and N is the
 * number of sequences in the database. The values beta, alpha, lambda and
 * K are statistical, Karlin-Altschul parameters.
 *
 * The value of the length adjustment computed by this routine, A,
 * will always be an integer smaller than the fixed point of
 * f(ell). Usually, it will be the largest such integer.  However, the
 * computed length adjustment, A, will also be so small that
 *
 *    K * (m - A) * (n - N * A) > min(m,n).
 *
 * Moreover, an iterative method is used to compute A, and under
 * unusual circumstances the iterative method may not converge.
 *
 * @param K      the statistical parameter K
 * @param logK   the natural logarithm of K
 * @param alpha_d_lambda    the ratio of the statistical parameters
 *                          alpha and lambda (for ungapped alignments, the
 *                          value 1/H should be used)
 * @param beta              the statistical parameter beta (for ungapped
 *                          alignments, beta == 0)
 * @param query_length      the length of the query sequence
 * @param db_length         the length of the database
 * @param db_num_seq        the number of sequences in the database
 * @param length_adjustment the computed value of the length adjustment [out]
 *
 * @return   0 if length_adjustment is known to be the largest integer less
 *           than the fixed point of f(ell); 1 otherwise.
 */

int32 BlastComputeLengthAdjustment(
  double K, double alpha_d_lambda, double beta,
  int32 query_length, int64 db_length, int32 db_num_seqs,
  int32& length_adjustment)
{
  double logK = log(K);

    int32 i;                     /* iteration index */
    const int32 maxits = 20;     /* maximum allowed iterations */
    double m = query_length;
    double n = static_cast<double>(db_length);
    double N = db_num_seqs;

    double ell;            /* A float value of the length adjustment */
    double ss;             /* effective size of the search space */
    double ell_min = 0, ell_max;   /* At each iteration i,
                                         * ell_min <= ell <= ell_max. */
    bool converged    = false;       /* True if the iteration converged */
    double ell_next = 0;   /* Value the variable ell takes at iteration
                                 * i + 1 */
    /* Choose ell_max to be the largest nonnegative value that satisfies
     *
     *    K * (m - ell) * (n - N * ell) > max(m,n)
     *
     * Use quadratic formula: 2 c /( - b + sqrt( b*b - 4 * a * c )) */
    { /* scope of a, mb, and c, the coefficients in the quadratic formula
       * (the variable mb is -b) */
        double a  = N;
        double mb = m * N + n;
        double c  = n * m - std::max(m, n) / K;

        if(c < 0) {
            length_adjustment = 0;
            return 1;
        } else {
            ell_max = 2 * c / (mb + sqrt(mb * mb - 4 * a * c));
        }
    } /* end scope of a, mb and c */

    for(i = 1; i <= maxits; i++) {      /* for all iteration indices */
        double ell_bar;    /* proposed next value of ell */
        ell      = ell_next;
        ss       = (m - ell) * (n - N * ell);
        ell_bar  = alpha_d_lambda * (logK + log(ss)) + beta;
        if(ell_bar >= ell) { /* ell is no bigger than the true fixed point */
            ell_min = ell;
            if(ell_bar - ell_min <= 1.0) {
                converged = true;
                break;
            }
            if(ell_min == ell_max) { /* There are no more points to check */
                break;
            }
        } else { /* else ell is greater than the true fixed point */
            ell_max = ell;
        }
        if(ell_min <= ell_bar && ell_bar <= ell_max) {
          /* ell_bar is in range. Accept it */
            ell_next = ell_bar;
        } else { /* else ell_bar is not in range. Reject it */
            ell_next = (i == 1) ? ell_max : (ell_min + ell_max) / 2;
        }
    } /* end for all iteration indices */
    if(converged) { /* the iteration converged */
        /* If ell_fixed is the (unknown) true fixed point, then we
         * wish to set length_adjustment to floor(ell_fixed).  We
         * assume that floor(ell_min) = floor(ell_fixed) */
        length_adjustment = (int32) ell_min;
        /* But verify that ceil(ell_min) != floor(ell_fixed) */
        ell = ceil(ell_min);
        if( ell <= ell_max ) {
          ss = (m - ell) * (n - N * ell);
          if(alpha_d_lambda * (logK + log(ss)) + beta >= ell) {
            /* ceil(ell_min) == floor(ell_fixed) */
            length_adjustment = (int32) ell;
          }
        }
    } else { /* else the iteration did not converge. */
        /* Use the best value seen so far */
        length_adjustment = (int32) ell_min;
    }

    return converged ? 0 : 1;
}

int32 BlastComputeLengthAdjustment(const Matrix& inMatrix, int32 query_length,
                                   int64 db_length, int32 db_num_seqs)
{
  int32 lengthAdjustment;
  (void)BlastComputeLengthAdjustment(
      inMatrix.GappedKappa(),
      inMatrix.GappedAlpha() / inMatrix.GappedLambda(),
      inMatrix.GappedBeta(),
      query_length, db_length, db_num_seqs, lengthAdjustment);
  return lengthAdjustment;
}

} // namespace ncbi

// --------------------------------------------------------------------

template<int WORDSIZE>
struct Word
{
  static const uint32 kMaxWordIndex;
  static const uint32 kMaxIndex;

  Word()
  {
    for (uint32 i = 0; i <= WORDSIZE; ++i)
      aa[i] = 0;
  }

  Word(const uint8* inSequence)
  {
    for (uint32 i = 0; i < WORDSIZE; ++i)
      aa[i] = inSequence[i];
    aa[WORDSIZE] = 0;
  }

  uint8& operator[](uint32 ix) { return aa[ix]; }
  const uint8* c_str() const { return aa; }
  size_t length() const { return WORDSIZE; }

  class PermutationIterator
  {
    public:
      PermutationIterator(Word inWord, const Matrix& inMatrix,
                          int32 inThreshold)
        : mWord(inWord), mIndex(0), mMatrix(inMatrix), mThreshold(inThreshold)
      {}

    bool Next(uint32& outIndex);

    private:
      Word mWord;
      uint32 mIndex;
      const Matrix& mMatrix;
      int32 mThreshold;
  };

  uint8 aa[WORDSIZE + 1];
};

template<>
const uint32 Word<2>::kMaxWordIndex = 0x0003FF;
template<>
const uint32 Word<2>::kMaxIndex = kAACount * kAACount;
template<>
const uint32 Word<3>::kMaxWordIndex = 0x007FFF;
template<>
const uint32 Word<3>::kMaxIndex = kAACount * kAACount * kAACount;
template<>
const uint32 Word<4>::kMaxWordIndex = 0x0FFFFF;
template<>
const uint32 Word<4>::kMaxIndex = kAACount * kAACount * kAACount * kAACount;

template<int WORDSIZE>
bool Word<WORDSIZE>::PermutationIterator::Next(uint32& outIndex)
{
  bool result = false;
  Word w;

  while (mIndex < kMaxIndex)
  {
    uint32 ix = mIndex;
    ++mIndex;

    int32 score = 0;
    outIndex = 0;

    for (uint32 i = 0; i < WORDSIZE; ++i)
    {
      uint32 resNr = ix % kAACount;
      w[i] = resNr;
      ix /= kAACount;
      score += mMatrix(mWord[i], w[i]);
      outIndex = outIndex << kBits | resNr;
    }

    if (score >= mThreshold)
    {
      result = true;
      break;
    }
  }

  return result;
}

template<int WORDSIZE>
class WordHitIterator
{
  static const uint32 kMask;

  struct Entry
  {
    uint16 mCount;
    uint16 mDataOffset;
  };

  public:

  typedef Word<WORDSIZE> IWord;
  typedef typename IWord::PermutationIterator WordPermutationIterator;

  struct WordHitIteratorStaticData
  {
    std::vector<Entry> mLookup;
    std::vector<uint16> mOffsets;
  };

  WordHitIterator(const WordHitIteratorStaticData& inStaticData)
  : mLookup(inStaticData.mLookup),
    mOffsets(inStaticData.mOffsets),
    mTargetCurrent(nullptr),
    mTargetEnd(nullptr),
    mIndex(0),
    mOffset(0),
    mCount(0),
    mTargetOffset(0)
  {
  }

  static void Init(const sequence& inQuery, const Matrix& inMatrix,
                   uint32 inThreshhold,
                   WordHitIteratorStaticData& outStaticData);

  void Reset(const sequence& inTarget);
  bool Next(uint16& outQueryOffset, uint16& outTargetOffset);
  uint32 Index() const { return mIndex; }

  private:

  const uint8* mTargetCurrent;
  const uint8* mTargetEnd;
  uint16 mTargetOffset;
  const std::vector<Entry>& mLookup;
  const std::vector<uint16>& mOffsets;
  uint32 mIndex;
  const uint16* mOffset;
  uint16 mCount;
};

template<> const uint32 WordHitIterator<2>::kMask = 0x0001F;
template<> const uint32 WordHitIterator<3>::kMask = 0x003FF;
template<> const uint32 WordHitIterator<4>::kMask = 0x07FFF;

template<int WORDSIZE>
void WordHitIterator<WORDSIZE>::Init(const sequence& inQuery,
                                     const Matrix& inMatrix,
                                     uint32 inThreshhold,
                                     WordHitIteratorStaticData& outStaticData)
{
  uint64 N = IWord::kMaxWordIndex;
  size_t M = 0;

  std::vector<std::vector<uint16>> test(N);

  for (uint16 i = 0; i < inQuery.length() - WORDSIZE + 1; ++i)
  {
    IWord w(inQuery.c_str() + i);

    WordPermutationIterator p(w, inMatrix, inThreshhold);
    uint32 ix;

    while (p.Next(ix))
    {
      test[ix].push_back(i);
      ++M;
    }
  }

  outStaticData.mLookup = std::vector<Entry>(N);
  outStaticData.mOffsets = std::vector<uint16>(M);

  uint16* data = &outStaticData.mOffsets[0];

  for (uint32 i = 0; i < N; ++i)
  {
    outStaticData.mLookup[i].mCount = static_cast<uint16>(test[i].size());
    outStaticData.mLookup[i].mDataOffset = static_cast<uint16>(
        data - &outStaticData.mOffsets[0]);

    for (uint32 j = 0; j < outStaticData.mLookup[i].mCount; ++j)
      *data++ = test[i][j];
  }

  assert(data == &outStaticData.mOffsets[0] + M);
#if not defined(NDEBUG)
  outStaticData.mOffsets.push_back(0);
#endif
}

template<int WORDSIZE>
void WordHitIterator<WORDSIZE>::Reset(const sequence& inTarget)
{
  mTargetCurrent = inTarget.c_str();
  mTargetEnd = mTargetCurrent + inTarget.length();
  mTargetOffset = 0;
  mIndex = 0;

  for (uint32 i = 0; i < WORDSIZE and mTargetCurrent != mTargetEnd; ++i)
    mIndex = mIndex << kBits | *mTargetCurrent++;

  Entry current = mLookup[mIndex];
  mCount = current.mCount;
  mOffset = &mOffsets[current.mDataOffset];
}

template<int WORDSIZE>
bool WordHitIterator<WORDSIZE>::Next(uint16& outQueryOffset,
                                     uint16& outTargetOffset)
{
  bool result = false;

  for (;;)
  {
    if (mCount-- > 0)
    {
      outQueryOffset = *mOffset++;
      outTargetOffset = mTargetOffset;
      result = true;
      break;
    }

    if (mTargetCurrent == mTargetEnd)
      break;

    mIndex = ((mIndex & kMask) << kBits) | *mTargetCurrent++;
    ++mTargetOffset;

    Entry current = mLookup[mIndex];
    mCount = current.mCount;
    mOffset = &mOffsets[current.mDataOffset];
  }

  return result;
}

// --------------------------------------------------------------------

struct DiagonalStartTable
{
  DiagonalStartTable()
  : mTable(nullptr),
    mTableLength(0),
    mTargetLength(0)
  {
  }

  ~DiagonalStartTable() { delete[] mTable; }

  void  Reset(int32 inQueryLength, int32 inTargetLength)
  {
    mTargetLength = inTargetLength;

    int32 n = inQueryLength + inTargetLength + 1;
    if (mTable == nullptr or n >= mTableLength)
    {
      uint32 k = ((n / 10240) + 1) * 10240;
      int32* t = new int32[k];
      delete[] mTable;
      mTable = t;
      mTableLength = k;
    }

    std::fill(mTable, mTable + n, -inTargetLength);
  }

  int32&  operator()(uint16 inQueryOffset, uint16 inTargetOffset)
  {
    return mTable[mTargetLength - inTargetOffset + inQueryOffset];
  }

  private:
    DiagonalStartTable(const DiagonalStartTable&);
    DiagonalStartTable&  operator=(const DiagonalStartTable&);

    int32* mTable;
    int32 mTableLength;
    int32 mTargetLength;
};

// --------------------------------------------------------------------

struct DPData
{
  DPData(size_t inDimX, size_t inDimY) : mDimX(inDimX), mDimY(inDimY)
  {
    mDPDataLength = (inDimX + 1) * (inDimY + 1);
    mDPData = new int16[mDPDataLength];
  }
  ~DPData() { delete[] mDPData; }

  int16    operator()(uint32 inI, uint32 inJ) const
  {
    return mDPData[inI * mDimY + inJ];
  }

  int16&    operator()(uint32 inI, uint32 inJ)
  {
    return mDPData[inI * mDimY + inJ];
  }

  int16* mDPData;
  size_t mDPDataLength;
  size_t mDimX;
  size_t mDimY;
};

struct DiscardTraceBack
{
  int16 operator()(int16 inB, int16 inIx, int16 inIy, uint32 /*inI*/,
                   uint32 /*inJ*/) const
  {
    return std::max(std::max(inB, inIx), inIy);
  }

  void Set(uint32 inI, uint32 inJ, int16 inD) {}
};

struct RecordTraceBack
{
  RecordTraceBack(DPData& inTraceBack) : mTraceBack(inTraceBack) {}

  int16 operator()(int16 inB, int16 inIx, int16 inIy, uint32 inI, uint32 inJ)
  {
    int16 result;

    if (inB >= inIx and inB >= inIy)
    {
      result = inB;
      mTraceBack(inI, inJ) = 0;
    }
    else if (inIx >= inB and inIx >= inIy)
    {
      result = inIx;
      mTraceBack(inI, inJ) = 1;
    }
    else
    {
      result = inIy;
      mTraceBack(inI, inJ) = -1;
    }

    return result;
  }

  void Set(uint32 inI, uint32 inJ, int16 inD) { mTraceBack(inI, inJ) = inD; }

  DPData&  mTraceBack;
};

// --------------------------------------------------------------------

inline void ReadEntry(const char*& inFasta, const char* inEnd,
                      sequence& outTarget)
{
  assert(inFasta == inEnd or *inFasta == '>');

  while (inFasta != inEnd and *inFasta++ != '\n');

  outTarget.clear();

  bool bol = false;
  while (inFasta != inEnd)
  {
    char ch = *inFasta++;

    if (ch == '\n')
      bol = true;
    else if (ch == '>' and bol)
    {
      --inFasta;
      break;
    }
    else
    {
      uint8 rn = ResidueNr(ch);
      if (rn < kResCount)
        outTarget += rn;
      bol = false;
    }
  }
}

// --------------------------------------------------------------------

struct Hsp
{
  uint32 mScore;
  uint32 mQueryStart;
  uint32 mQueryEnd;
  uint32 mTargetStart;
  uint32 mTargetEnd;
  sequence mAlignedQuery;
  sequence mAlignedTarget;
  double mBitScore;
  double mExpect;
  bool mGapped;

  bool operator>(const Hsp& inHsp) const
  {
    return mScore > inHsp.mScore;
  }

  void CalculateExpect(int64 inSearchSpace, double inLambda,
                       double inLogKappa);
  bool Overlaps(const Hsp& inOther) const
  {
    return mQueryEnd >= inOther.mQueryStart and
           mQueryStart <= inOther.mQueryEnd and
           mTargetEnd >= inOther.mTargetStart and
           mTargetStart <= inOther.mTargetEnd;
  }
};

void Hsp::CalculateExpect(int64 inSearchSpace, double inLambda,
                          double inLogKappa)
{
  mBitScore = floor((inLambda * mScore - inLogKappa) / kLn2);
  mExpect = inSearchSpace / pow(2., mBitScore);
}

// --------------------------------------------------------------------

struct Hit;
typedef std::shared_ptr<Hit> HitPtr;

struct Hit
{
  Hit(const char* inEntry, const sequence& inTarget);

  void AddHsp(const Hsp& inHsp);
  void Cleanup(int64 inSearchSpace, double inLambda, double inLogKappa,
               double inExpect);

  std::string mDefLine;
  sequence mTarget;
  std::vector<Hsp> mHsps;
};

Hit::Hit(const char* inEntry, const sequence& inTarget)
  : mDefLine(inEntry, const_cast<const char*>(strchr(inEntry, '\n'))),
    mTarget(inTarget)
{
}

void Hit::AddHsp(const Hsp& inHsp)
{
  bool found = false;

  foreach (auto& hsp, mHsps)
  {
    if (inHsp.Overlaps(hsp))
    {
      if (hsp.mScore < inHsp.mScore)
        hsp = inHsp;
      found = true;
      break;
    }
  }

  if (not found)
    mHsps.push_back(inHsp);
}

void Hit::Cleanup(int64 inSearchSpace, double inLambda, double inLogKappa,
                  double inExpect)
{
  std::sort(mHsps.begin(), mHsps.end(), std::greater<Hsp>());

  std::vector<Hsp>::iterator a = mHsps.begin();
  while (a != mHsps.end() and a + 1 != mHsps.end())
  {
    std::vector<Hsp>::iterator b = a + 1;
    while (b != mHsps.end())
    {
      if (a->Overlaps(*b))
        b = mHsps.erase(b);
      else
        ++b;
    }
    ++a;
  }

  for_each(mHsps.begin(), mHsps.end(), [=](Hsp& hsp) {
    hsp.CalculateExpect(inSearchSpace, inLambda, inLogKappa);
  });

  std::sort(mHsps.begin(), mHsps.end(), std::greater<Hsp>());

  mHsps.erase(
    remove_if(mHsps.begin(), mHsps.end(), [=](const Hsp& hsp) -> bool {
      return hsp.mExpect > inExpect;
    }),
    mHsps.end());
}

// --------------------------------------------------------------------

template<int WORDSIZE>
class BlastQuery
{
  public:
    BlastQuery(const std::string& inQuery, bool inFilter, double inExpect,
               const std::string& inMatrix, bool inGapped, int32 inGapOpen,
               int32 inGapExtend, uint32 inReportLimit);
    ~BlastQuery();

    void Search(const std::vector<fs::path>& inDatabanks,
                MProgress& inProgress, uint32 inNrOfThreads);
    void WriteAsFasta(std::ostream& inStream);

  private:
    void SearchPart(const char* inFasta, size_t inLength,
                    MProgress& inProgress, uint32& outDbCount,
                    int64& outDbLength, std::vector<HitPtr>& outHits) const;

    int32 Extend(int32& ioQueryStart, const sequence& inTarget,
                      int32& ioTargetStart, int32& ioDistance) const;
    template<class Iterator1, class Iterator2, class TraceBack>
    int32 AlignGapped(Iterator1 inQueryBegin, Iterator1 inQueryEnd,
                      Iterator2 inTargetBegin, Iterator2 inTargetEnd,
                      TraceBack& inTraceBack, int32 inDropOff,
                      uint32& outBestX, uint32& outBestY) const;

    int32 AlignGappedFirst(const sequence& inTarget, Hsp& ioHsp) const;
    int32 AlignGappedSecond(const sequence& inTarget, Hsp& ioHsp) const;

    void AddHit(HitPtr inHit, std::vector<HitPtr>& inHitList) const;

    typedef WordHitIterator<WORDSIZE> IWordHitIterator;
    typedef typename IWordHitIterator::WordHitIteratorStaticData StaticData;

    std::string mUnfiltered;
    sequence mQuery;
    Matrix mMatrix;
    double mExpect;
    double mCutOff;
    bool mGapped;
    int32 mS1;
    int32 mS2;
    int32 mXu;
    int32 mXg;
    int32 mXgFinal;
    uint32 mReportLimit;

    uint32 mDbCount;
    int64 mDbLength;
    int64 mSearchSpace;

    std::vector<HitPtr> mHits;

    StaticData mWordHitData;
};

template<int WORDSIZE>
BlastQuery<WORDSIZE>::BlastQuery(const std::string& inQuery, bool inFilter,
                                 double inExpect, const std::string& inMatrix,
                                 bool inGapped, int32 inGapOpen,
                                 int32 inGapExtend, uint32 inReportLimit)
  : mUnfiltered(inQuery),
    mMatrix(inMatrix, inGapOpen, inGapExtend),
    mExpect(inExpect),
    mGapped(inGapped),
    mReportLimit(inReportLimit),
    mDbCount(0),
    mDbLength(0),
    mSearchSpace(0),
    mCutOff(0)
{
  if (mQuery.length() >= kMaxSequenceLength)
    throw mas_exception("Query length exceeds maximum");

  mUnfiltered.erase(remove_if(mUnfiltered.begin(),
                              mUnfiltered.end(), [](char aa) -> bool {
    return ResidueNr(aa) >= kResCount;
  }), mUnfiltered.end());

  std::string query(mUnfiltered);
  if (inFilter)
    query = filter::SEG(query);

  transform(query.begin(), query.end(),
            back_inserter(mQuery), [](char aa) -> uint8 {
    return ResidueNr(aa);
  });

  mXu = static_cast<int32>(
      ceil((kLn2 * kUngappedDropOff) / mMatrix.UngappedLambda()));
  mXg = static_cast<int32>(
      (kLn2 * kGappedDropOff) / mMatrix.GappedLambda());
  mXgFinal = static_cast<int32>(
      (kLn2 * kGappedDropOffFinal) / mMatrix.GappedLambda());
  mS1 = static_cast<int32>(
      (kLn2 * kGapTrigger + log(mMatrix.UngappedKappa())) / mMatrix.UngappedLambda());

  // we're not using S2
  mS2 = static_cast<int32>(
      (kLn2 * kGapTrigger + log(mMatrix.GappedKappa())) / mMatrix.GappedLambda());

  IWordHitIterator::Init(mQuery, mMatrix, kThreshold, mWordHitData);
}

template<int WORDSIZE>
BlastQuery<WORDSIZE>::~BlastQuery()
{
}

template<int WORDSIZE>
void BlastQuery<WORDSIZE>::Search(const std::vector<fs::path>& inDatabanks,
                                  MProgress& inProgress, uint32 inNrOfThreads)
{
  foreach (const fs::path& p, inDatabanks)
  {
    io::mapped_file file(p.string().c_str(), io::mapped_file::readonly);
    if (not file.is_open())
      throw mas_exception(boost::format("FastA file %s not open") % p);

    const char* data = file.const_data();
    size_t length = file.size();

    if (inNrOfThreads <= 1)
      SearchPart(data, length, inProgress, mDbCount, mDbLength, mHits);
    else
    {
      boost::thread_group t;
      boost::mutex m;

      size_t k = length / inNrOfThreads;
      for (uint32 i = 0; i < inNrOfThreads and length > 0; ++i)
      {
        size_t n = k;
        if (n > length)
          n = length;
        const char* end = data + n;
        while (n < length and *end != '>')
          ++end, ++n;

        t.create_thread([data, n, &m, &inProgress, this]() {
          uint32 dbCount = 0;
          int64 dbLength = 0;
          std::vector<HitPtr> hits;

          this->SearchPart(data, n, inProgress, dbCount, dbLength, hits);

          boost::mutex::scoped_lock lock(m);
          mDbCount += dbCount;
          mDbLength += dbLength;
          this->mHits.insert(mHits.end(), hits.begin(), hits.end());
        });

        data += n;
        length -= n;
      }

      t.join_all();
    }
  }

  int32 lengthAdjustment = ncbi::BlastComputeLengthAdjustment(
      mMatrix, static_cast<uint32>(mQuery.length()), mDbLength, mDbCount);

  int64 effectiveQueryLength = mQuery.length() - lengthAdjustment;
  int64 effectiveDbLength = mDbLength - mDbCount * lengthAdjustment;

  mSearchSpace = effectiveDbLength * effectiveQueryLength;

  if (not mHits.empty())
  {
    boost::thread_group t;
    boost::detail::atomic_count ix(-1);

    for (uint32 i = 0; i < inNrOfThreads; ++i)
    {
      t.create_thread([this, &ix]() {
        double lambda = mMatrix.GappedLambda();
        double logK = log(mMatrix.GappedKappa());

        for (;;)
        {
          uint32 next = ++ix;
          if (next >= mHits.size())
            break;

          HitPtr hit = mHits[next];

          foreach (Hsp& hsp, hit->mHsps)
            hsp.mScore = this->AlignGappedSecond(hit->mTarget, hsp);

          hit->Cleanup(mSearchSpace, lambda, logK, mExpect);
        }
      });
    }
    t.join_all();
  }

  mHits.erase(
    remove_if(mHits.begin(), mHits.end(),
              [](const HitPtr hit) -> bool { return hit->mHsps.empty(); }),
    mHits.end());

  std::sort(mHits.begin(), mHits.end(),
            [](const HitPtr a, const HitPtr b) -> bool {
    return a->mHsps.front().mScore > b->mHsps.front().mScore or
           (a->mHsps.front().mScore == b->mHsps.front().mScore and
            a->mDefLine < b->mDefLine);
  });

  if (mHits.size() > mReportLimit and mReportLimit > 0)
    mHits.erase(mHits.begin() + mReportLimit, mHits.end());
}


template<int WORDSIZE>
void BlastQuery<WORDSIZE>::WriteAsFasta(std::ostream& inStream)
{
  foreach (HitPtr hit, mHits)
  {
    std::string seq;
    foreach (uint8 r, hit->mTarget)
    {
      if (seq.length() % 73 == 72)
        seq += '\n';
      seq += kResidues[r];
    }

    inStream << hit->mDefLine << std::endl
         << seq << std::endl;
  }
}

template<int WORDSIZE>
void BlastQuery<WORDSIZE>::SearchPart(const char* inFasta, size_t inLength,
                                      MProgress& inProgress,
  uint32& outDbCount, int64& outDbLength, std::vector<HitPtr>& outHits) const
{
  const char* end = inFasta + inLength;
  int32 queryLength = static_cast<int32>(mQuery.length());

  IWordHitIterator iter(mWordHitData);
  DiagonalStartTable diagonals;
  sequence target;
  target.reserve(kMaxSequenceLength);

  int64 hitsToDb = 0, extensions = 0, successfulExtensions = 0;
  HitPtr hit;

  while (inFasta != end)
  {
    if (hit)
    {
      AddHit(hit, outHits);
      hit.reset();
    }

    const char* entry = inFasta;
    ReadEntry(inFasta, end, target);

    inProgress.Consumed(inFasta - entry);

    if (target.empty() or target.length() > kMaxSequenceLength)
      continue;

    outDbCount += 1;
    outDbLength += target.length();

    iter.Reset(target);
    diagonals.Reset(queryLength, static_cast<int32>(target.length()));

    uint16 queryOffset, targetOffset;
    while (iter.Next(queryOffset, targetOffset))
    {
      ++hitsToDb;

      int32& ds = diagonals(queryOffset, targetOffset);
      int32 distance = queryOffset - ds;

      if (distance >= kHitWindow)
        ds = queryOffset;
      else if (distance > WORDSIZE)
      {
        int32 queryStart = ds;
        int32 targetStart = targetOffset - distance;
        int32 alignmentDistance = distance + WORDSIZE;

        if (targetStart < 0 or queryStart < 0)
          continue;

        ++extensions;

        int32 score = Extend(queryStart, target, targetStart,
                             alignmentDistance);

        if (score >= mS1)
        {
          ++successfulExtensions;

          Hsp hsp;

          // extension results, to be updated later
          hsp.mQueryStart = queryStart;
          hsp.mQueryEnd = queryStart + alignmentDistance;
          hsp.mTargetStart = targetStart;
          hsp.mTargetEnd = targetStart + alignmentDistance;

          if (not hit)
            hit.reset(new Hit(entry, target));

          if (mGapped)
            hsp.mScore = AlignGappedFirst(target, hsp);
          else
          {
            hsp.mScore = score;
            hsp.mAlignedQuery = mQuery.substr(
                hsp.mQueryStart, hsp.mQueryEnd - hsp.mQueryStart);
            hsp.mAlignedTarget = hit->mTarget.substr(
                hsp.mTargetStart, hsp.mTargetEnd - hsp.mTargetStart);
          }

          hit->AddHsp(hsp);
        }

        ds = queryStart + alignmentDistance;
      }
    }
  }

  if (hit)
    AddHit(hit, outHits);
}

template<int WORDSIZE>
int32 BlastQuery<WORDSIZE>::Extend(int32& ioQueryStart,
                                   const sequence& inTarget,
                                   int32& ioTargetStart,
                                   int32& ioDistance) const
{
  // use iterators
  sequence::const_iterator ai = mQuery.begin() + ioQueryStart;
  sequence::const_iterator bi = inTarget.begin() + ioTargetStart;

  int32 score = 0;
  for (int i = 0; i < ioDistance; ++i, ++ai, ++bi)
    score += mMatrix(*ai, *bi);

  // record start and stop positions for optimal score
  sequence::const_iterator qe = ai;

  for (int32 test = score,
       n = static_cast<int32>(
         std::min(mQuery.end() - ai, inTarget.end() - bi));
       test >= score - mXu and n > 0;
       --n, ++ai, ++bi)
  {
    test += mMatrix(*ai, *bi);

    if (test > score)
    {
      score = test;
      qe = ai;
    }
  }

  ai = mQuery.begin() + ioQueryStart;
  bi = inTarget.begin() + ioTargetStart;
  sequence::const_iterator qs = ai + 1;

  for (int32 test = score, n = std::min(ioQueryStart, ioTargetStart);
     test >= score - mXu and n > 0;
     --n)
  {
    test += mMatrix(*--ai, *--bi);

    if (test > score)
    {
      score = test;
      qs = ai;
    }
  }

  int32 delta = static_cast<int32>(ioQueryStart - (qs - mQuery.begin()));
  ioQueryStart -= delta;
  ioTargetStart -= delta;
  ioDistance = static_cast<int32>(qe - qs);

  return score;
}

template<int WORDSIZE>
template<class Iterator1, class Iterator2, class TraceBack>
int32 BlastQuery<WORDSIZE>::AlignGapped(
    Iterator1 inQueryBegin, Iterator1 inQueryEnd, Iterator2 inTargetBegin,
    Iterator2 inTargetEnd, TraceBack& inTraceBack, int32 inDropOff,
    uint32& outBestX, uint32& outBestY) const
{
  const Matrix& s = mMatrix;  // for readability
  TraceBack& tb_max = inTraceBack;
  int32 d = s.OpenCost();
  int32 e = s.ExtendCost();

  uint32 dimX = static_cast<uint32>(inQueryEnd - inQueryBegin);
  uint32 dimY = static_cast<uint32>(inTargetEnd - inTargetBegin);

  DPData B(dimX, dimY);
  DPData Ix(dimX, dimY);
  DPData Iy(dimX, dimY);

  int32 bestScore = 0;
  uint32 bestX;
  uint32 bestY;
  uint32 colStart = 1;
  uint32 lastColStart = 1;
  uint32 colEnd = dimY;

  // first column
  uint32 i = 1, j = 1;
  Iterator1 x = inQueryBegin;
  Iterator2 y = inTargetBegin;

  // first cell
  int32 Ix1 = kSentinalScore, Iy1 = kSentinalScore;

  // (1)
  int32 M = s(*x, *y);

  // (2)
  (void)tb_max(M, kSentinalScore, kSentinalScore, 1, 1);
  bestScore = B(1, 1) = M;
  bestX = bestY = 1;

  // (3)
  Ix(1, 1) = M - d;

  // (4)
  Iy(1, 1) = M - d;

  // remaining cells in the first column
  y = inTargetBegin + 1;
  Ix1 = kSentinalScore;
  M = kSentinalScore;

  for (j = 2; y != inTargetEnd; ++j, ++y)
  {
    Iy1 = Iy(i, j - 1);

    // (2)
    int32 Bij = B(i, j) = Iy1;
    tb_max.Set(i, j, -1);

    // (3)
    Ix(i, j) = kSentinalScore;

    // (4)
    Iy(i, j) = Iy1 - e;

    if (Bij < bestScore - inDropOff)
    {
      colEnd = j;
      break;
    }
  }

  // remaining columns
  ++x;
  for (i = 2; x != inQueryEnd and colEnd >= colStart; ++i, ++x)
  {
    y = inTargetBegin + colStart - 1;
    uint32 newColStart = colStart;
    bool beforeFirstRow = true;

    for (j = colStart; y != inTargetEnd; ++j, ++y)
    {
      Ix1 = kSentinalScore;
      Iy1 = kSentinalScore;

      if (j < colEnd)
        Ix1 = Ix(i - 1, j);

      if (j > colStart)
        Iy1 = Iy(i, j - 1);

      // (1)
      if (j <= lastColStart or j > colEnd)
        M = kSentinalScore;
      else
        M = B(i - 1, j - 1) + s(*x, *y);

      // cut off the max value
      if (M > std::numeric_limits<int16>::max())
        M = std::numeric_limits<int16>::max();

      // (2)
      int32 Bij = B(i, j) = tb_max(M, Ix1, Iy1, i, j);

      // (3)
      Ix(i, j) = std::max(M - d, Ix1 - e);

      // (4)
      Iy(i, j) = std::max(M - d, Iy1 - e);

      if (Bij > bestScore)
      {
        bestScore = Bij;
        bestX = i;
        bestY = j;
        beforeFirstRow = false;
      }
      else if (Bij < bestScore - inDropOff)
      {
        if (beforeFirstRow)
        {
          newColStart = j;
          if (newColStart > colEnd)
            break;
        }
        else if (j > bestY + 1)
        {
          colEnd = j;
          break;
        }
      }
      else
      {
        beforeFirstRow = false;
        if (j > colEnd)
          colEnd = j;
      }
    }

    lastColStart = colStart;
    colStart = newColStart;
  }

  outBestY = bestY;
  outBestX = bestX;

  return bestScore;
}

template<int WORDSIZE>
int32 BlastQuery<WORDSIZE>::AlignGappedFirst(const sequence& inTarget,
                                             Hsp& ioHsp) const
{
  int32 score;

  uint32 x, y;
  uint32 targetSeed = (ioHsp.mTargetStart + ioHsp.mTargetEnd) / 2;
  uint32 querySeed = (ioHsp.mQueryStart + ioHsp.mQueryEnd) / 2;

  DiscardTraceBack tb;

  score = AlignGapped(
    mQuery.begin() + querySeed + 1, mQuery.end(),
    inTarget.begin() + targetSeed + 1, inTarget.end(),
    tb, mXg, x, y);

  score += AlignGapped(
    mQuery.rbegin() + (mQuery.length() - querySeed), mQuery.rend(),
    inTarget.rbegin() + (inTarget.length() - targetSeed), inTarget.rend(),
    tb, mXg, x, y);

  score += mMatrix(mQuery[querySeed], inTarget[targetSeed]);

  return score;
}

template<int WORDSIZE>
int32 BlastQuery<WORDSIZE>::AlignGappedSecond(const sequence& inTarget,
                                              Hsp& ioHsp) const
{
  uint32 x, y;
  int32 score = 0;

  ioHsp.mGapped = false;

  sequence alignedQuery;
  sequence alignedTarget;

  uint32 targetSeed = (ioHsp.mTargetStart + ioHsp.mTargetEnd) / 2;
  uint32 querySeed = (ioHsp.mQueryStart + ioHsp.mQueryEnd) / 2;

  // start with the part before the seed
  DPData d1(querySeed + 1, targetSeed + 1);
  RecordTraceBack tbb(d1);

  score = AlignGapped(
    mQuery.rbegin() + (mQuery.length() - querySeed), mQuery.rend(),
    inTarget.rbegin() + (inTarget.length() - targetSeed), inTarget.rend(),
    tbb, mXgFinal, x, y);
  ioHsp.mQueryStart = querySeed - x;
  ioHsp.mTargetStart = targetSeed - y;

  sequence::const_iterator qi = mQuery.begin() + querySeed - x, qis = qi;
  sequence::const_iterator si = inTarget.begin() + targetSeed - y, sis = si;

  uint32 qLen = 1;
  uint32 sLen = 1;

  while (x >= 1 and y >= 1)
  {
    if (x >= 1 and y >= 1 and d1(x, y) == 0)
    {
      alignedQuery += *qi++;
      alignedTarget += *si++;
      --x;
      --y;
    }
    else if (y >= 1 and d1(x, y) < 0)
    {
      alignedQuery += '-';
      alignedTarget += *si++;
      --y;
      ioHsp.mGapped = true;
    }
    else // if (x >= 1 and d1(x, y) > 0)
    {
      alignedQuery += *qi++;
      alignedTarget += '-';
      --x;
      ioHsp.mGapped = true;
    }
  }

  qLen += static_cast<uint32>(qi - qis);
  sLen += static_cast<uint32>(si - sis);

  // the seed itself
  alignedQuery += mQuery[querySeed];
  alignedTarget += inTarget[targetSeed];
  score += mMatrix(mQuery[querySeed], inTarget[targetSeed]);

  // and the part after the seed
  DPData d2(mQuery.length() - querySeed, inTarget.length() - targetSeed);
  RecordTraceBack tba(d2);

  score += AlignGapped(
    mQuery.begin() + querySeed + 1, mQuery.end(),
    inTarget.begin() + targetSeed + 1, inTarget.end(),
    tba, mXgFinal, x, y);

  sequence::const_reverse_iterator qri = mQuery.rbegin() + (mQuery.length() - querySeed) - 1 - x, qris = qri;
  sequence::const_reverse_iterator sri = inTarget.rbegin() + (inTarget.length() - targetSeed) - 1 - y, sris = sri;

  sequence q, s;

  while (x >= 1 and y >= 1)
  {
    if (x >= 1 and y >= 1 and d2(x, y) == 0)
    {
      q += *qri++;
      s += *sri++;
      --x;
      --y;
    }
    else if (y >= 1 and d2(x, y) < 0)
    {
      q += '-';
      s += *sri++;
      --y;
      ioHsp.mGapped = true;
    }
    else // if (x >= 1 and d2(x, y) > 0)
    {
      q += *qri++;
      s += '-';
      --x;
      ioHsp.mGapped = true;
    }
  }

  reverse(q.begin(), q.end());
  reverse(s.begin(), s.end());

  alignedQuery += q;
  alignedTarget += s;

  qLen += static_cast<uint32>(qri - qris);
  sLen += static_cast<uint32>(sri - sris);

  ioHsp.mAlignedQuery.assign(alignedQuery.begin(), alignedQuery.end());
  ioHsp.mAlignedTarget.assign(alignedTarget.begin(), alignedTarget.end());
  ioHsp.mQueryEnd = ioHsp.mQueryStart + qLen;
  ioHsp.mTargetEnd = ioHsp.mTargetStart + sLen;

  return score;
}

template<int WORDSIZE>
void BlastQuery<WORDSIZE>::AddHit(HitPtr inHit,
                                  std::vector<HitPtr>& inHitList) const
{
  std::sort(inHit->mHsps.begin(), inHit->mHsps.end(), std::greater<Hsp>());

  inHitList.push_back(inHit);

  auto cmp = [](const HitPtr a, const HitPtr b) -> bool {
    return a->mHsps.front().mScore > b->mHsps.front().mScore;
  };

  push_heap(inHitList.begin(), inHitList.end(), cmp);
  if (inHitList.size() > mReportLimit and mReportLimit > 0)
  {
    pop_heap(inHitList.begin(), inHitList.end(), cmp);
    inHitList.erase(inHitList.end() - 1);
  }
}


void SearchAndWriteResultsAsFastA(
    std::ostream& inOutFile, const std::vector<fs::path>& inDatabanks,
    const std::string& inQuery, const std::string& inProgram,
    const std::string& inMatrix, uint32 inWordSize, double inExpect,
    bool inFilter, bool inGapped, int32 inGapOpen, int32 inGapExtend,
    uint32 inReportLimit, uint32 inThreads)
{
  if (inProgram != "blastp")
    throw mas_exception(boost::format("Unsupported program %s") % inProgram);

  if (inGapped)
  {
    if (inGapOpen == -1) inGapOpen = 11;
    if (inGapExtend == -1) inGapExtend = 1;
  }

  if (inWordSize == 0) inWordSize = 3;

  std::string query(inQuery), queryID("query"), queryDef;

  if (ba::starts_with(inQuery, ">"))
  {
    boost::smatch m;
    if (regex_search(inQuery, m, kFastARE, boost::match_not_dot_newline))
    {
      queryID = m[4];
      if (queryID.empty())
        queryID = m[2];
      queryDef = m[7];
      query = m.suffix();
    }
    else
    {
      queryID = inQuery.substr(1, inQuery.find('\n') - 1);
      query = inQuery.substr(queryID.length() + 2, std::string::npos);

      std::string::size_type s = queryID.find(' ');
      if (s != std::string::npos)
      {
        queryDef = queryID.substr(s + 1);
        queryID.erase(s, std::string::npos);
      }
    }
  }

  int64 totalLength = accumulate(inDatabanks.begin(), inDatabanks.end(), 0LL,
    [](int64 l, const fs::path& p) -> int64 { return l + fs::file_size(p); });

  MProgress progress(totalLength, "blast");

  switch (inWordSize)
  {
    case 2:
    {
      BlastQuery<2> q(query, inFilter, inExpect, inMatrix, inGapped, inGapOpen,
                      inGapExtend, inReportLimit);
      q.Search(inDatabanks, progress, inThreads);
      q.WriteAsFasta(inOutFile);
      break;
    }
    case 3:
    {
      BlastQuery<3> q(query, inFilter, inExpect, inMatrix, inGapped, inGapOpen,
                      inGapExtend, inReportLimit);
      q.Search(inDatabanks, progress, inThreads);
      q.WriteAsFasta(inOutFile);
      break;
    }
    case 4:
    {
      BlastQuery<4> q(query, inFilter, inExpect, inMatrix, inGapped, inGapOpen,
                      inGapExtend, inReportLimit);
      q.Search(inDatabanks, progress, inThreads);
      q.WriteAsFasta(inOutFile);
      break;
    }
    default:
      throw mas_exception(
          boost::format("Unsupported word size %d") % inWordSize);
  }
}
