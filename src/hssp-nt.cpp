#include "mas.h"

#include <iostream>
#include <set>
#include <cmath>

#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/filesystem/operations.hpp>
#include <boost/filesystem/convenience.hpp>
#include <boost/filesystem/fstream.hpp>
#include <boost/format.hpp>
#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/filter/zlib.hpp>
#include <boost/iostreams/device/back_inserter.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/copy.hpp>
#include <boost/iostreams/filter/bzip2.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/foreach.hpp>
#define foreach BOOST_FOREACH
#include <boost/program_options.hpp>
#include <boost/program_options/config.hpp>
#include <boost/date_time/gregorian/gregorian.hpp>
#include <boost/date_time/date_clock_device.hpp>
#include <boost/regex.hpp>
#include <boost/tr1/cmath.hpp>

#include "utils.h"
#include "structure.h"
#include "dssp.h"
#include "matrix.h"
#include "buffer.h"

using namespace std;
namespace fs = boost::filesystem;
namespace ba = boost::algorithm;
namespace po = boost::program_options;
namespace io = boost::iostreams;

int VERBOSE = 0;

// --------------------------------------------------------------------

namespace HSSP
{

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

bool drop(uint32 ident, uint32 length, float threshold)
{
	uint32 ix = max(10U, min(length, 80U)) - 10;
	return ident < length * (kHomologyThreshold[ix] + threshold);
}

// --------------------------------------------------------------------

const int8 kMPam250[] = {
	  2,                                                                                                               // A
	  0,   3,                                                                                                          // B
	 -2,  -4,  12,                                                                                                     // C
	  0,   3,  -5,   4,                                                                                                // D
	  0,   3,  -5,   3,   4,                                                                                           // E
	 -3,  -4,  -4,  -6,  -5,   9,                                                                                      // F
	  1,   0,  -3,   1,   0,  -5,   5,                                                                                 // G
	 -1,   1,  -3,   1,   1,  -2,  -2,   6,                                                                            // H
	 -1,  -2,  -2,  -2,  -2,   1,  -3,  -2,   5,                                                                       // I
	 -1,   1,  -5,   0,   0,  -5,  -2,   0,  -2,   5,                                                                  // K
	 -2,  -3,  -6,  -4,  -3,   2,  -4,  -2,   2,  -3,   6,                                                             // L
	 -1,  -2,  -5,  -3,  -2,   0,  -3,  -2,   2,   0,   4,   6,                                                        // M
	  0,   2,  -4,   2,   1,  -3,   0,   2,  -2,   1,  -3,  -2,   2,                                                   // N
	  1,  -1,  -3,  -1,  -1,  -5,   0,   0,  -2,  -1,  -3,  -2,   0,   6,                                              // P
	  0,   1,  -5,   2,   2,  -5,  -1,   3,  -2,   1,  -2,  -1,   1,   0,   4,                                         // Q
	 -2,  -1,  -4,  -1,  -1,  -4,  -3,   2,  -2,   3,  -3,   0,   0,   0,   1,   6,                                    // R
	  1,   0,   0,   0,   0,  -3,   1,  -1,  -1,   0,  -3,  -2,   1,   1,  -1,   0,   2,                               // S
	  1,   0,  -2,   0,   0,  -3,   0,  -1,   0,   0,  -2,  -1,   0,   0,  -1,  -1,   1,   3,                          // T
	  0,  -2,  -2,  -2,  -2,  -1,  -1,  -2,   4,  -2,   2,   2,  -2,  -1,  -2,  -2,  -1,   0,   4,                     // V
	 -6,  -5,  -8,  -7,  -7,   0,  -7,  -3,  -5,  -3,  -2,  -4,  -4,  -6,  -5,   2,  -2,  -5,  -6,  17,                // W
	 -3,  -3,   0,  -4,  -4,   7,  -5,   0,  -1,  -4,  -1,  -2,  -2,  -5,  -4,  -4,  -3,  -3,  -2,   0,  10,           // Y
	  0,   2,  -5,   3,   3,  -5,   0,   2,  -2,   0,  -3,  -2,   1,   0,   3,   0,   0,  -1,  -2,  -6,  -4,   3,      // Z
	  0,  -1,  -3,  -1,  -1,  -2,  -1,  -1,  -1,  -1,  -1,  -1,   0,  -1,  -1,  -1,   0,   0,  -1,  -4,  -2,  -1,  -1, // x
};

const float
	kMPam250ScalingFactor = log(2.f) / 3.0f,
	kMPam250MisMatchAverage = -1.484210526f;

// Dayhoff matrix as used by maxhom
const float kDayhoffData[] =
{
	 1.5f,																																 // A
	 0.0f, 0.0f,                                                                                                                         // B
	 0.3f, 0.0f, 1.5f,                                                                                                                   // C
	 0.3f, 0.0f,-0.5f, 1.5f,                                                                                                             // D
	 0.3f, 0.0f,-0.6f, 1.0f, 1.5f,                                                                                                       // E
	-0.5f, 0.0f,-0.1f,-1.0f,-0.7f, 1.5f,                                                                                                 // F
	 0.7f, 0.0f, 0.2f, 0.7f, 0.5f,-0.6f, 1.5f,                                                                                           // G
	-0.1f, 0.0f,-0.1f, 0.4f, 0.4f,-0.1f,-0.2f, 1.5f,                                                                                     // H
	 0.0f, 0.0f, 0.2f,-0.2f,-0.2f, 0.7f,-0.3f,-0.3f, 1.5f,                                                                               // I
	 0.0f, 0.0f,-0.6f, 0.3f, 0.3f,-0.7f,-0.1f, 0.1f,-0.2f, 1.5f,                                                                         // K
	-0.1f, 0.0f,-0.8f,-0.5f,-0.3f, 1.2f,-0.5f,-0.2f, 0.8f,-0.3f, 1.5f,                                                                   // L
	 0.0f, 0.0f,-0.6f,-0.4f,-0.2f, 0.5f,-0.3f,-0.3f, 0.6f, 0.2f, 1.3f, 1.5f,                                                             // M
	 0.2f, 0.0f,-0.3f, 0.7f, 0.5f,-0.5f, 0.4f, 0.5f,-0.3f, 0.4f,-0.4f,-0.3f, 1.5f,                                                       // N
	 0.5f, 0.0f, 0.1f, 0.1f, 0.1f,-0.7f, 0.3f, 0.2f,-0.2f, 0.1f,-0.3f,-0.2f, 0.0f, 1.5f,                                                 // P
	 0.2f, 0.0f,-0.6f, 0.7f, 0.7f,-0.8f, 0.2f, 0.7f,-0.3f, 0.4f,-0.1f, 0.0f, 0.4f, 0.3f, 1.5f,                                           // Q
	-0.3f, 0.0f,-0.3f, 0.0f, 0.0f,-0.5f,-0.3f, 0.5f,-0.3f, 0.8f,-0.4f, 0.2f, 0.1f, 0.3f, 0.4f, 1.5f,                                     // R
	 0.4f, 0.0f, 0.7f, 0.2f, 0.2f,-0.3f, 0.6f,-0.2f,-0.1f, 0.2f,-0.4f,-0.3f, 0.3f, 0.4f,-0.1f, 0.1f, 1.5f,                               // S
	 0.4f, 0.0f, 0.2f, 0.2f, 0.2f,-0.3f, 0.4f,-0.1f, 0.2f, 0.2f,-0.1f, 0.0f, 0.2f, 0.3f,-0.1f,-0.1f, 0.3f, 1.5f,                         // T
	 0.2f, 0.0f, 0.2f,-0.2f,-0.2f, 0.2f, 0.2f,-0.3f, 1.1f,-0.2f, 0.8f, 0.6f,-0.3f, 0.1f,-0.2f,-0.3f,-0.1f, 0.2f, 1.5f,                   // V
	-0.8f, 0.0f,-1.2f,-1.1f,-1.1f, 1.3f,-1.0f,-0.1f,-0.5f, 0.1f, 0.5f,-0.3f,-0.3f,-0.8f,-0.5f, 1.4f, 0.3f,-0.6f,-0.8f, 1.5f,             // W
	-0.3f, 0.0f, 1.0f,-0.5f,-0.5f, 1.4f,-0.7f, 0.3f, 0.1f,-0.6f, 0.3f,-0.1f,-0.1f,-0.8f,-0.6f,-0.6f,-0.4f,-0.3f,-0.1f, 1.1f, 1.5f,       // Y
	 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f  // Z
};

// 22 real letters and 1 dummy
const char kResidues[] = "ABCDEFGHIKLMNPQRSTVWYZX";

inline uint8 ResidueNr(char inAA)
{
	int result = 23;

	const static uint8 kResidueNrTable[] = {
	//	A   B   C   D   E   F   G   H   I       K   L   M   N       P   Q   R   S   T  U=X  V   W   X   Y   Z
		0,  1,  2,  3,  4,  5,  6,  7,  8, 23,  9, 10, 11, 12, 23, 13, 14, 15, 16, 17, 22, 18, 19, 22, 20, 21
	};
	
	inAA |= 040;
	if (inAA >= 'a' and inAA <= 'z')
		result = kResidueNrTable[inAA - 'a'];
	return result;
}

template<typename T>
inline T score(const T inMatrix[], uint8 inAA1, uint8 inAA2)
{
	T result;

	if (inAA1 >= inAA2)
		result = inMatrix[(inAA1 * (inAA1 + 1)) / 2 + inAA2];
	else
		result = inMatrix[(inAA2 * (inAA2 + 1)) / 2 + inAA1];

	return result;	
}

sequence encode(const string& s)
{
	sequence result(s.length(), 0);
	transform(s.begin(), s.end(), result.begin(), [](char aa) -> uint8 { return aa == '-' ? '-' : ResidueNr(aa); });
	return result;
}

string decode(const sequence& s)
{
	string result(s.length(), 0);
	transform(s.begin(), s.end(), result.begin(), [](uint8 r) -> char { return r == '-' ? '-' : kResidues[r]; });
	return result;
}

bool is_gap(char aa)
{
	return aa == ' ' or aa == '.' or aa == '-';
}

bool is_hydrophilic(char aa)
{
	aa &= ~040;
	return aa == 'D' or aa == 'E' or aa == 'G' or aa == 'K' or aa == 'N' or aa == 'Q' or aa == 'P' or aa == 'R' or aa == 'S';
}

// --------------------------------------------------------------------

float calculateDistance(const sequence& a, const sequence& b)
{
	const float kDistanceGapOpen = 10;
	const float kDistanceGapExtend = 0.2f;

	int32 x = 0, dimX = static_cast<int32>(a.length());
	int32 y = 0, dimY = static_cast<int32>(b.length());

	matrix<float>	B(dimX, dimY);
	matrix<float>	Ix(dimX, dimY);
	matrix<float>	Iy(dimX, dimY);
	matrix<uint16>	id(dimX, dimY);
	
	Ix(0, 0) = 0;
	Iy(0, 0) = 0;
	
	uint16 highId = 0;
	
	Ix(x, y) = 0;
	Iy(x, y) = 0;
	if (x > 0 and y > 0)
		id(x - 1, y - 1) = highId;

	int32 startX = x, startY = y;
	float high = -numeric_limits<float>::max();
	uint16 highIdSub = 0;

	for (x = startX; x < dimX; ++x)
	{
		for (y = startY; y < dimY; ++y)
		{
			float Ix1 = 0; if (x > startX) Ix1 = Ix(x - 1, y);
			float Iy1 = 0; if (y > startY) Iy1 = Iy(x, y - 1);

			// (1)
			float M = score(kMPam250, a[x], b[y]);
			if (x > startX and y > startY)
				M += B(x - 1, y - 1);

			float s;
			uint16 i = 0;
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
			Ix(x, y) = max(M - kDistanceGapOpen, Ix1 - kDistanceGapExtend);
			
			// (4)
			Iy(x, y) = max(M - kDistanceGapOpen, Iy1 - kDistanceGapExtend);
		}
	}

	highId += highIdSub;
	
	float result = 1.0f - float(highId) / max(dimX, dimY);

	assert(result >= 0.0f);
	assert(result <= 1.0f);
	
	return result;
}

// --------------------------------------------------------------------

struct MHit
{
	MHit(const MHit& e);
	
	MHit(const string& id, const string& def, const sequence& seq)
		: m_id(id), m_def(def), m_seq(seq), m_distance(0)
		, m_identical(0), m_similar(0), m_length(0), m_gaps(0), m_gapn(0) {}

	struct insertion
	{
		uint32			m_ipos, m_jpos;
		string			m_seq;
	};
	
	static MHit*		Create(const string& id, const string& def,
							const string& seq, const sequence& chain);

	void				Update(const sequence& inChain);
	void				Update(const matrix<int8>& inTraceBack, const sequence& inChain,
							int32 inX, int32 inY, const matrix<float>& inB);

	string				m_id, m_acc, m_def;
	sequence			m_seq;
	string				m_aligned;
	float				m_distance, m_score;
	int32				m_ifir, m_ilas, m_jfir, m_jlas;
	uint32				m_identical, m_similar, m_length;
	uint32				m_gaps, m_gapn;
	vector<insertion>	m_insertions;
};

ostream& operator<<(ostream& os, const MHit& hit)
{
	string seq = hit.m_aligned;
	foreach (char& r, seq)
		if (r == ' ' or r == '.') r = '-';
	
	for (string::size_type i = 72; i < seq.length(); i += 73)
		seq.insert(seq.begin() + i, '\n');
	
	os << '>' << hit.m_id /*<< ' ' << hit.m_def*/ << endl
	   << seq << endl;
	
	return os;
}

MHit* MHit::Create(const string& id, const string& def, const string& seq, const sequence& chain)
{
	MHit* result = new MHit(id, def, encode(seq));

	static const boost::regex
		kM6FastARE("^(\\w+)((?:\\|([^| ]*))(?:\\|([^| ]+))?(?:\\|([^| ]+))?(?:\\|([^| ]+))?)");
	
	boost::smatch m;
	if (boost::regex_match(result->m_id, m, kM6FastARE, boost::match_not_dot_newline))
	{
		if (m[1] == "sp" or m[1] == "tr")
		{
			result->m_acc = m[3];
			result->m_id = m[4];
		}
		else
			result->m_id = m[2];
	}
	
	result->m_aligned = string(chain.length(), ' ');
	result->m_distance = calculateDistance(result->m_seq, chain);
	return result;
}

void MHit::Update(const sequence& inChain)
{
	assert(inChain.length() == m_aligned.length());
	
	bool gap = false;

	uint32 x = 0;
	while (is_gap(m_aligned[x]))
		++x;
	uint32 lx = x;
	
	insertion ins;
	ins.m_ipos = m_ifir;
	ins.m_jpos = m_jfir;
	
	for ( ; x < m_aligned.length(); ++x)
	{
		bool igap = is_gap(inChain[x]);
		bool jgap = is_gap(m_aligned[x]);
		
		if (not igap)
			++ins.m_ipos;

		if (not jgap)
			++ins.m_jpos;

		if (igap and jgap)
			continue;
		
		if (igap)
		{
			if (not gap)
			{
				m_aligned[lx] |= 040;
				m_insertions.push_back(ins);
				m_insertions.back().m_seq += m_aligned[lx];
			}
			gap = true;
			m_insertions.back().m_seq += m_aligned[x];
		}
		else if (not jgap)
		{
			lx = x;
			if (gap)
			{
				m_aligned[x] |= 040;
				m_insertions.back().m_seq += m_aligned[x];
				gap = false;
			}
		}
	}
}

void MHit::Update(const matrix<int8>& inTraceBack, const sequence& inChain,
	int32 inX, int32 inY, const matrix<float>& inB)
{
	m_ilas = inX + 1;
	m_jlas = inY + 1;

	int32 x = inX;
	int32 y = inY;
	
	bool gap = false;
	
	insertion ins;
	
	// trace back the matrix
	while (x >= 0 and y >= 0 and inB(x, y) > 0)
	{
		++m_length;
		switch (inTraceBack(x, y))
		{
			case -1:
				if (not gap)
				{
					++m_gaps;
					int32 nx = x + 1;
					while (nx < m_aligned.length() and is_gap(m_aligned[nx]))
						++nx;
					assert(nx < m_aligned.length());
					m_aligned[nx] |= 040;
					
					ins.m_seq = m_aligned[nx];
				}					
				++m_gapn;
				ins.m_seq += kResidues[m_seq[y]];
				gap = true;
				--y;
				break;

			case 1:
				m_aligned[x] = '.';
				--x;
				break;

			case 0:
				m_aligned[x] = kResidues[m_seq[y]];
				if (gap)
				{
					m_aligned[x] |= 040;
					ins.m_seq += m_aligned[x];
					ins.m_ipos = x + 1;
					ins.m_jpos = y + 1;
					
					reverse(ins.m_seq.begin(), ins.m_seq.end());
					m_insertions.push_back(ins);
				}
				gap = false;

				if (inChain[x] == m_seq[y])
					++m_identical, ++m_similar;
				else if (score(kMPam250, inChain[x], m_seq[y]) > 0)
					++m_similar;
				--x;
				--y;
				break;

			default:
				assert(false);
		}
	}
	
	assert(gap == false);
	m_ifir = x + 2;
	m_jfir = y + 2;
	
	m_score = float(m_identical) / m_length;
}

// --------------------------------------------------------------------

struct MResInfo
{
	uint8			m_letter;
	uint8			m_chain_id;
	uint32			m_seq_nr;
	uint32			m_pdb_nr;
	MSecondaryStructure
					m_ss;
	string			m_dssp;
	float			m_consweight;
	uint32			m_nocc;
	uint32			m_dist[23];
	float			m_dist_weight[23];
	float			m_sum_dist_weight;
	uint32			m_ins, m_del;
	float			m_score[23];

	void			Add(uint8 r, float inDistance);
};

typedef vector<MResInfo> MResInfoList;

void MResInfo::Add(uint8 r, float inDistance)
{
	float weight = 1 - inDistance;
	
	assert(r < 23);
	m_nocc += 1;
	m_dist[r] += 1;
	
	m_dist_weight[r] += weight;
	m_sum_dist_weight += weight;
	
	for (int i = 0; i < 23; ++i)
	{
		float si = 0;
		
		for (int j = 0; j < 23; ++j)
			si += float(score(kMPam250, i, j)) * m_dist_weight[j];
		
		m_score[i] = si / m_sum_dist_weight;
	}
}

// --------------------------------------------------------------------

struct MProfile
{
					MProfile(const MChain& inChain, const sequence& inSequence,
						float inThreshold, float inFragmentCutOff);

	void			Process(istream& inHits, progress& inProgress);
	void			Align(MHit* e);
	void			AdjustGapCosts(vector<float>& gop, vector<float>& gep);

	void			dump(ostream& os, const matrix<int8>& tb, const sequence& s);

	void			PrintFastA();

	const MChain&	m_chain;
	sequence		m_seq;
	MResInfoList	m_residues;
	vector<MHit*>	m_entries;
	float			m_threshold, m_frag_cutoff;
};

MProfile::MProfile(const MChain& inChain, const sequence& inSequence, float inThreshold, float inFragmentCutOff)
	: m_chain(inChain), m_seq(inSequence)
	, m_threshold(inThreshold), m_frag_cutoff(inFragmentCutOff)
{
	const vector<MResidue*>& residues = m_chain.GetResidues();
	vector<MResidue*>::const_iterator ri = residues.begin();

	uint32 seq_nr = 1;
	for (uint32 i = 0; i < inSequence.length(); ++i)
	{
		assert(ri != residues.end());
		
		if (ri != residues.begin() and (*ri)->GetNumber() > (*(ri - 1))->GetNumber() + 1)
			++seq_nr;

		string dssp = ResidueToDSSPLine(**ri).substr(5, 34);
		MResInfo res = { inSequence[i], m_chain.GetChainID(), seq_nr,
			(*ri)->GetNumber(), (*ri)->GetSecondaryStructure(), dssp };
		res.Add(res.m_letter, 0);
		m_residues.push_back(res);
		
		++ri;
		++seq_nr;
	}
}

void MProfile::dump(ostream& os, const matrix<int8>& tb, const sequence& s)
{
	os << ' ';
	foreach (auto& e, m_residues)
		os << kResidues[e.m_letter];
	os << endl;

	for (uint32 y = 0; y < s.length(); ++y)
	{
		os << kResidues[s[y]];
		for (uint32 x = 0; x < m_residues.size(); ++x)
		{
			switch (tb(x, y))
			{
				case -1:	os << '|'; break;
				case 0:		os << '\\'; break;
				case 1:		os << '-'; break;
				case 2:		os << '.'; break;
			}
		}
		os << endl;
	}
}

const float kResidueSpecificPenalty[22] = {
	1.13f,		// A
	1.00f,		// B
	1.13f,		// C
	0.96f,		// D
	1.31f,		// E
	1.20f,		// F
	0.61f,		// G
	1.00f,		// H
	1.32f,		// I
	0.96f,		// K
	1.21f,		// L
	1.29f,		// M
	0.63f,		// N
	0.74f,		// P
	1.07f,		// Q
	0.72f,		// R
	0.76f,		// S
	0.89f,		// T
	1.25f,		// V
	1.23f,		// W
	1.00f,		// Y
	1.00f,		// Z
};

void MProfile::AdjustGapCosts(vector<float>& gop, vector<float>& gep)
{
	assert(gop.size() == m_seq.length());
	assert(gop.size() == m_residues.size());

	for (int32 ix = 0; ix < m_residues.size(); ++ix)
	{
		MResInfo& e = m_residues[ix];
		float resSpecific = 1.0f;
		
		switch (e.m_ss)
		{
			case alphahelix:
			case helix_5:
			case helix_3:
				resSpecific = 1.25;
				break;

			case betabridge:
			case strand:
				resSpecific = 1.25;
				break;

			default:
				resSpecific = kResidueSpecificPenalty[e.m_letter];
				break;
		}

		// if there is a gap, lower gap open cost
		if (e.m_del > 0)
		{
			gop[ix] *= 0.3f * ((1.0f + m_entries.size() - e.m_del) / (m_entries.size() + 1));
			gep[ix] /= 2;
		}
		
		// else if there is a gap within 8 residues, increase gap cost
		else
		{
			for (int32 d = 0; d < 8; ++d)
			{
				if (ix + d >= int32(m_residues.size()) or
					(m_residues[ix + d].m_del + m_residues[ix + d].m_ins) > 0 or
					ix - d < 0 or
					(m_residues[ix - d].m_del + m_residues[ix - d].m_ins) > 0)
				{
					gop[ix] *= (2 + ((8 - d) * 2)) / 8.f;
					break;
				}
			}
			
			gop[ix] *= resSpecific;
		}
	}
}

void MProfile::Align(MHit* e)
{
	int32 x = 0, dimX = static_cast<int32>(m_seq.length());
	int32 y = 0, dimY = static_cast<int32>(e->m_seq.length());
	
	matrix<float> B(dimX, dimY);
	matrix<float> Ix(dimX, dimY);
	matrix<float> Iy(dimX, dimY);
	matrix<int8> tb(dimX, dimY);
	
	float minLength = static_cast<float>(dimX), maxLength = static_cast<float>(dimY);
	if (minLength > maxLength)
		swap(minLength, maxLength);
	
	float gop = 10, gep = 0.2f;
	
	float logmin = 1.0f / log10(minLength);
	float logdiff = 1.0f + 0.5f * log10(minLength / maxLength);
	
	// initial gap open cost, 0.05f is the remaining magical number here...
	float magic = 1; //0.05f;
	gop = (gop / (logdiff * logmin)) * abs(kMPam250MisMatchAverage) * kMPam250ScalingFactor * magic;

	// position specific gap penalties
	// initial gap extend cost is adjusted for difference in sequence lengths
	vector<float> gop_a(dimX, gop), gep_a(dimX, gep /* * (1 + log10(float(dimX) / dimY))*/);
	AdjustGapCosts(gop_a, gep_a);
	
	vector<float> gop_b(dimY, gop), gep_b(dimY, gep /* * (1 + log10(float(dimY) / dimX))*/);
//	adjust_gp(gop_b, gep_b, b);

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

//cerr << x << ',' << y << ' '
//	 << kResidues[m_residues[x].m_letter] << kResidues[e->m_seq[y]] << ' '
//	 << Ix1 << ',' << Iy1 << ',' << M << " => ";

			float s;
			if (M >= Ix1 and M >= Iy1)
			{
				tb(x, y) = 0;
				B(x, y) = s = M;
				
				Ix(x, y) = M - (x < dimX - 1 ? gop_a[x] : 0);
				Iy(x, y) = M - (y < dimY - 1 ? gop_b[y] : 0);
			}
			else if (Ix1 >= Iy1)
			{
				tb(x, y) = 1;
				B(x, y) = s = Ix1;

				Ix(x, y) = Ix1 - gep_a[x];
				Iy(x, y) = max(M - (y < dimY - 1 ? gop_b[y] : 0), Iy1 - gep_b[y]);
			}
			else
			{
				tb(x, y) = -1;
				B(x, y) = s = Iy1;

				Ix(x, y) = max(M - (x < dimX - 1 ? gop_a[x] : 0), Ix1 - gep_a[x]);
				Iy(x, y) = Iy1 - gep_b[y];
			}
			
			if (highS < s)
			{
				highS = s;
				highX = x;
				highY = y;
			}

//cerr << Ix(x, y) << ',' << Iy(x, y) << ',' << B(x, y) << " => " << int(tb(x, y)) << endl;
		}
	}

//#ifndef NDEBUG
//
//ofstream log("alignment.log");
//if (not log.is_open()) throw mas_exception("open log");
//
//log << "B" << endl;
//for (int y = 0; y < dimY; ++y)
//{
//	for (int x = 0; x < dimX; ++x)
//	{
//		log << B(x, y) << '\t';
//	}
//	log << endl;
//}
//	
//log << "Ix" << endl;
//for (int y = 0; y < dimY; ++y)
//{
//	for (int x = 0; x < dimX; ++x)
//	{
//		log << Ix(x, y) << '\t';
//	}
//	log << endl;
//}
//	
//log << "Iy" << endl;
//for (int y = 0; y < dimY; ++y)
//{
//	for (int x = 0; x < dimX; ++x)
//	{
//		log << Iy(x, y) << '\t';
//	}
//	log << endl;
//}
//	
//log << "tb" << endl;
//for (int y = 0; y < dimY; ++y)
//{
//	for (int x = 0; x < dimX; ++x)
//	{
//		log << int(tb(x, y)) << '\t';
//	}
//	log << endl;
//}
//
//#endif


	// build the alignment
	x = highX;
	y = highY;

//	if (VERBOSE >= 6)
//		dump(cerr, tb, e->m_seq);

	uint32 ident = 0, length = 0, xgaps = 0;

	// trace back the matrix
	while (x >= 0 and y >= 0 and B(x, y) > 0)
	{
		++length;
		switch (tb(x, y))
		{
			case -1:
				--y;
				++xgaps;
				break;

			case 1:
				--x;
				break;

			case 0:
				if (e->m_seq[y] == m_seq[x])
				{
					if (m_seq[x] == '-')
						--length;
					else
						++ident;
				}
				--x;
				--y;
				break;
			
			default:
				assert(false);
				break;
		}
	}
	
	if (m_seq.length() * m_frag_cutoff < length and not drop(ident, length, m_threshold))
	{
#if CUT
		e->m_distance = 1 - float(ident) / length;
		e->Update(tb, m_seq, highX, highY, B);
		m_entries.push_back(e);
		for (uint32 i = 0; i < m_seq.length(); ++i)
		{
			if (is_gap(e->m_aligned[i]))
				continue;
			m_residues[i].Add(ResidueNr(e->m_aligned[i]), e->m_distance);
		}
		
		// update insert/delete counters for the residues
		x = highX;
		y = highY;
		bool gap = false;
		while (x >= 0 and y >= 0 and B(x, y) > 0)
		{
			switch (tb(x, y))
			{
				case -1:
					if (not gap)
						++m_residues[x + 1].m_ins;
					gap = true;
					--y;
					break;
	
				case 1:
					++m_residues[x].m_del;
					--x;
					break;
	
				case 0:
					gap = false;
					--x;
					--y;
					break;
			}
		}
#else
		// Add the hit since it is within the required parameters.
		// Calculate the new distance
		e->m_distance = 1 - float(ident) / length;

		// reserve space, if needed
		if (xgaps > 0)
		{
			uint32 n = (((m_residues.size() + xgaps) / 256) + 1) * 256;
			m_seq.reserve(n);
			foreach (MHit* e, m_entries)
				e->m_aligned.reserve(n);
		}
		
		// update insert/delete counters for the residues
		x = highX;	e->m_ilas = x + 1;
		y = highY;	e->m_jlas = y + 1;
		bool gap = false;
		e->m_aligned = string(m_seq.length() + xgaps, ' ');
		
		const static MResInfo rgap = {};
		
		while (x >= 0 and y >= 0 and B(x, y) > 0)
		{
			++e->m_length;
			
			switch (tb(x, y))
			{
				case -1:
					e->m_aligned[x + xgaps] = kResidues[e->m_seq[y]];

					if (not gap)
					{
						++m_residues[x].m_ins;
						++e->m_gaps;
					}					
					++e->m_gapn;
					gap = true;
					
					m_residues.insert(m_residues.begin() + x + 1, rgap);
					m_residues[x + 1].m_del = m_entries.size() + 1;
//					m_residues[x + 1].Add(e->m_seq[y], e->m_distance);
					m_seq.insert(m_seq.begin() + x + 1, '-');
					
					foreach (MHit* e, m_entries)
						e->m_aligned.insert(e->m_aligned.begin() + x + 1, '-');
					
					--y;
					--xgaps;
					break;
	
				case 1:
					e->m_aligned[x + xgaps] = '-';
					++m_residues[x + 1].m_del;
					--x;
					break;
	
				case 0:
					e->m_aligned[x + xgaps] = kResidues[e->m_seq[y]];
					m_residues[x].Add(e->m_seq[y], e->m_distance);

					if (gap)
						e->m_aligned[x + xgaps] |= 040;
					gap = false;

					if (m_seq[x] == e->m_seq[y])
						++e->m_identical, ++e->m_similar;
					else if (score(kMPam250, m_seq[x], e->m_seq[y]) > 0)
						++e->m_similar;

					--x;
					--y;
					break;
			}
		}

		// update the new entry

		assert(gap == false);
		e->m_ifir = x + 2;
		e->m_jfir = y + 2;
		
		e->m_score = float(e->m_identical) / e->m_length;

		e->Update(m_seq);

		m_entries.push_back(e);
#endif

		//PrintFastA();
	}
	else
		delete e;
}

void MProfile::PrintFastA()
{
	string s = decode(m_seq);
	for (string::size_type i = 72; i < s.length(); i += 73)
		s.insert(s.begin() + i, '\n');
	
	cout << '>' << "PDB" << endl
		 << s << endl;
	
	foreach (MHit* e, m_entries)
		cout << *e;
}

// --------------------------------------------------------------------

void MProfile::Process(istream& inHits, progress& inProgress)
{
	string id, def, seq;
	for (;;)
	{
		string line;
		getline(inHits, line);
		if (line.empty() and inHits.eof())
			break;

		inProgress.step(line.length() + 1);
		
		if (ba::starts_with(line, ">"))
		{
			if (not (id.empty() or seq.empty()))
				Align(MHit::Create(id, def, seq, m_seq));
			
			id.clear();
			def.clear();
			seq.clear();
			
			string::size_type s = line.find(' ');
			if (s != string::npos)
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
		Align(MHit::Create(id, def, seq, m_seq));
	
	if (VERBOSE)
		PrintFastA();

	sort(m_entries.begin(), m_entries.end(), [](const MHit* a, const MHit* b) -> bool {
		return a->m_score > b->m_score;
	});
}

// --------------------------------------------------------------------

// Find the minimal set of overlapping sequences
// In case of strong similarity (distance <= 0.01) we take the longest chain.
void ClusterSequences(vector<sequence>& s, vector<uint32>& ix)
{
	for (;;)
	{
		bool found = false;
		for (uint32 i = 0; not found and i < s.size() - 1; ++i)
		{
			for (uint32 j = i + 1; not found and j < s.size(); ++j)
			{
				sequence& a = s[i];
				sequence& b = s[j];

				if (a.empty() or b.empty())
					continue;
				
				if (a == b)
				{
					s[j].clear();
					ix[j] = i;
					found = true;
				}
				else
				{
					float d = calculateDistance(a, b);
					// rescale distance to shortest length:
					d = 1 - (1 - d) * max(a.length(), b.length()) / min(a.length(), b.length());
					if (d <= 0.01)
					{
						if (b.length() > a.length())
							swap(i, j);

						s[j].clear();
						ix[j] = i;
						found = true;
					}
				}
			}
		}
		
		if (not found)
			break;
	}
}

// --------------------------------------------------------------------
// Write collected information as a HSSP file to the output stream

void CreateHSSPOutput(const string& inProteinID, const string& inProteinDescription,
	float inThreshold, float inFragmentCutOff, uint32 inSeqLength, uint32 inNChain, uint32 inKChain,
	const string& inUsedChains, const vector<MHit*>& inHits, const vector<MResInfo>& inResInfo, ostream& os)
{
	using namespace boost::gregorian;
	date today = day_clock::local_day();

	// print the header
	os << "HSSP       HOMOLOGY DERIVED SECONDARY STRUCTURE OF PROTEINS, VERSION 3.0 2012" << endl
	   << "PDBID      " << inProteinID << endl
	   << "DATE       file generated on " << to_iso_extended_string(today) << endl
//	   << "SEQBASE    " << inDatabank->GetName() << " version " << inDatabank->GetVersion() << endl
	   << "THRESHOLD  according to: t(L)=(290.15 * L ** -0.562) + " << (inThreshold * 100) << endl
	   << "CUTOFF     minimal alignment length as fraction of chain length: " << inFragmentCutOff << endl
	   << "REFERENCE  Sander C., Schneider R. : Database of homology-derived protein structures. Proteins, 9:56-68 (1991)." << endl
	   << "CONTACT    Maintained at http://www.cmbi.ru.nl/ by Maarten L. Hekkelman <m.hekkelman@cmbi.ru.nl>" << endl
	   << inProteinDescription
	   << boost::format("SEQLENGTH %5.5d") % inSeqLength << endl
	   << boost::format("NCHAIN     %4.4d chain(s) in %s data set") % inNChain % inProteinID << endl;
	
	if (inKChain != inNChain)
		os << boost::format("KCHAIN     %4.4d chain(s) used here ; chains(s) : ") % inKChain << inUsedChains << endl;
	
	os << boost::format("NALIGN     %4.4d") % inHits.size() << endl
	   << "NOTATION : ID: EMBL/SWISSPROT identifier of the aligned (homologous) protein" << endl
	   << "NOTATION : STRID: if the 3-D structure of the aligned protein is known, then STRID is the Protein Data Bank identifier as taken" << endl
	   << "NOTATION :   from the database reference or DR-line of the EMBL/SWISSPROT entry" << endl
	   << "NOTATION : %IDE: percentage of residue identity of the alignment" << endl
	   << "NOTATION : %SIM (%WSIM):  (weighted) similarity of the alignment" << endl
	   << "NOTATION : IFIR/ILAS: first and last residue of the alignment in the test sequence" << endl
	   << "NOTATION : JFIR/JLAS: first and last residue of the alignment in the alignend protein" << endl
	   << "NOTATION : LALI: length of the alignment excluding insertions and deletions" << endl
	   << "NOTATION : NGAP: number of insertions and deletions in the alignment" << endl
	   << "NOTATION : LGAP: total length of all insertions and deletions" << endl
	   << "NOTATION : LSEQ2: length of the entire sequence of the aligned protein" << endl
	   << "NOTATION : ACCNUM: SwissProt accession number" << endl
	   << "NOTATION : PROTEIN: one-line description of aligned protein" << endl
	   << "NOTATION : SeqNo,PDBNo,AA,STRUCTURE,BP1,BP2,ACC: sequential and PDB residue numbers, amino acid (lower case = Cys), secondary" << endl
	   << "NOTATION :   structure, bridge partners, solvent exposure as in DSSP (Kabsch and Sander, Biopolymers 22, 2577-2637(1983)" << endl
	   << "NOTATION : VAR: sequence variability on a scale of 0-100 as derived from the NALIGN alignments" << endl
	   << "NOTATION :   pair of lower case characters (AvaK) in the alignend sequence bracket a point of insertion in this sequence" << endl
	   << "NOTATION :   dots (....) in the alignend sequence indicate points of deletion in this sequence" << endl
	   << "NOTATION : SEQUENCE PROFILE: relative frequency of an amino acid type at each position. Asx and Glx are in their" << endl
	   << "NOTATION :   acid/amide form in proportion to their database frequencies" << endl
	   << "NOTATION : NOCC: number of aligned sequences spanning this position (including the test sequence)" << endl
	   << "NOTATION : NDEL: number of sequences with a deletion in the test protein at this position" << endl
	   << "NOTATION : NINS: number of sequences with an insertion in the test protein at this position" << endl
	   << "NOTATION : ENTROPY: entropy measure of sequence variability at this position" << endl
	   << "NOTATION : RELENT: relative entropy, i.e.  entropy normalized to the range 0-100" << endl
	   << "NOTATION : WEIGHT: conservation weight" << endl
	   << endl
	   << "## PROTEINS : identifier and alignment statistics" << endl
	   << "  NR.    ID         STRID   %IDE %WSIM IFIR ILAS JFIR JLAS LALI NGAP LGAP LSEQ2 ACCNUM     PROTEIN" << endl;
	   
	// print the first list
	uint32 nr = 1;
	boost::format fmt1("%5.5d : %12.12s%4.4s    %4.2f  %4.2f%5.5d%5.5d%5.5d%5.5d%5.5d%5.5d%5.5d%5.5d  %10.10s %s");
	foreach (MHit* h, inHits)
	{
		string id = h->m_id, acc = h->m_acc, pdb;

		if (id.length() > 12)
			id.erase(12, string::npos);
		else if (id.length() < 12)
			id.append(12 - id.length(), ' ');
		
		if (acc.length() > 10)
			acc.erase(10, string::npos);
		else if (acc.length() < 10)
			acc.append(10 - acc.length(), ' ');
		
		float sim = float(h->m_similar) / h->m_length;
		
		os << fmt1 % nr % id % pdb % h->m_score % sim
				   % h->m_ifir % h->m_ilas % h->m_jfir % h->m_jlas % h->m_length
				   % h->m_gaps % h->m_gapn % h->m_seq.length() % acc % h->m_def
		   << endl;
		
		++nr;
	}

	// print the alignments
	for (uint32 i = 0; i < inHits.size(); i += 70)
	{
		uint32 n = i + 70;
		if (n > inHits.size())
			n = inHits.size();
		
		uint32 k[7] = {
			((i +  0) / 10 + 1) % 10,
			((i + 10) / 10 + 1) % 10,
			((i + 20) / 10 + 1) % 10,
			((i + 30) / 10 + 1) % 10,
			((i + 40) / 10 + 1) % 10,
			((i + 50) / 10 + 1) % 10,
			((i + 60) / 10 + 1) % 10
		};
		
		os << boost::format("## ALIGNMENTS %4.4d - %4.4d") % (i + 1) % n << endl
		   << boost::format(" SeqNo  PDBNo AA STRUCTURE BP1 BP2  ACC NOCC  VAR  ....:....%1.1d....:....%1.1d....:....%1.1d....:....%1.1d....:....%1.1d....:....%1.1d....:....%1.1d")
		   					% k[0] % k[1] % k[2] % k[3] % k[4] % k[5] % k[6] << endl;

		int32 x = -1, nextNr = inResInfo.front().m_seq_nr;
		foreach (auto& ri, inResInfo)
		{
			++x;
			if (ri.m_chain_id == 0)
				continue;
				
			if (ri.m_seq_nr != nextNr)
				os << boost::format(" %5.5d        !  !           0   0    0    0    0") % ri.m_seq_nr << endl;

			string aln;
			
			foreach (MHit* hit, boost::make_iterator_range(inHits.begin() + i, inHits.begin() + n))
				aln += hit->m_aligned[x];
			
			uint32 ivar = uint32(100 * (1 - ri.m_consweight));

			os << ' ' << boost::format("%5.5d%s%4.4d %4.4d  ")
				% ri.m_seq_nr % ri.m_dssp % ri.m_nocc % ivar << aln << endl;

			nextNr = ri.m_seq_nr + 1;
		}
	}
	
	// ## SEQUENCE PROFILE AND ENTROPY
	os << "## SEQUENCE PROFILE AND ENTROPY" << endl
	   << " SeqNo PDBNo   V   L   I   M   F   W   Y   G   A   P   S   T   C   H   R   K   Q   E   N   D  NOCC NDEL NINS ENTROPY RELENT WEIGHT" << endl;
	
	int32 nextNr = inResInfo.front().m_seq_nr;
	foreach (auto& ri, inResInfo)
	{
		if (ri.m_chain_id == 0)
			continue;
		
		if (ri.m_seq_nr != nextNr)
			os << boost::format("%5.5d          0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0     0    0    0   0.000      0  1.00")
				% nextNr << endl;

		os << boost::format("%5.5d%5.5d %c") % ri.m_seq_nr % ri.m_pdb_nr % ri.m_chain_id;

		double entropy = 0;

		static const int8 kResIx[] = {
			//	V   L   I   M   F   W   Y   G   A   P   S   T   C   H   R   K   Q   E   N   D
			   18, 10,  8, 11,  5, 19, 20,  6,  0, 13, 16, 17,  2,  7, 15,  9, 14,  4, 12,  3
		};

		for (uint32 i = 0; i < 20; ++i)
		{
			double freq = double(ri.m_dist[kResIx[i]]) / ri.m_nocc;
			os << boost::format("%4.4d") % uint32(100.0 * freq + 0.5);
			if (freq > 0)
				entropy -= freq * log(freq);
		}

		uint32 relent = uint32(100 * entropy / log(20.0));
		os << "  " << boost::format("%4.4d %4.4d %4.4d   %5.3f   %4.4d  %4.2f") % ri.m_nocc % ri.m_del % ri.m_ins % entropy % relent % ri.m_consweight << endl;
		
		nextNr = ri.m_seq_nr + 1;
	}
	
	// insertion list
	
	os << "## INSERTION LIST" << endl
	   << " AliNo  IPOS  JPOS   Len Sequence" << endl;

	nr = 0;
	foreach (MHit* h, inHits)
	{
		++nr;
		foreach (auto& ins, h->m_insertions)
		{
			string s = ins.m_seq;
			
			uint32 ipos = inResInfo[ins.m_ipos].m_seq_nr;
			uint32 jpos = ins.m_jpos;

			os << boost::format(" %5.5d %5.5d %5.5d %5.5d ") % nr % ins.m_ipos % ins.m_jpos % (ins.m_seq.length() - 2);
			
			if (s.length() <= 100)
				os << s << endl;
			else
			{
				os << s.substr(0, 100) << endl;
				s.erase(0, 100);
				
				while (not s.empty())
				{
					uint32 n = s.length();
					if (n > 100)
						n = 100;
					
					os << "     +                   " << s.substr(0, n) << endl;
					s.erase(0, n);
				}
			}
		}			
	}
	
	os << "//" << endl;
}

// --------------------------------------------------------------------
// Calculate the variability of a residue, based on dayhoff similarity
// and weights

pair<const char*,uint32> kSentinel((const char*)nullptr, 0);

void CalculateConservation(buffer<pair<const char*,uint32>>& b,
	const vector<MHit*>& inHits, vector<float>& sumvar, vector<float>& sumdist)
{
	uint32 length = sumvar.size();
	vector<float> simval(length);

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
				simval[k] = numeric_limits<float>::min();

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
				if (simval[k] != numeric_limits<float>::min())
				{
					sumvar[k] += distance * simval[k];
					sumdist[k] += distance * 1.5f;
				}
			}
		}
	}

	b.put(kSentinel);
}

void CalculateConservation(const sequence& inChain, vector<MHit*>& inHits, MResInfoList& inResidues)
{
	if (VERBOSE)
		cerr << "Calculating conservation weights...";

	vector<float> sumvar(inChain.length(), 0), sumdist(inChain.length(), 0);
	
	// Calculate conservation weights in multiple threads to gain speed.
	buffer<pair<const char*,uint32>> b;
	boost::thread_group threads;
	boost::mutex sumLock;

	for (uint32 t = 0; t < boost::thread::hardware_concurrency(); ++t)
		threads.create_thread([&]() {
			vector<float> csumvar(sumvar.size(), 0), csumdist(sumdist.size(), 0);
			
			CalculateConservation(b, inHits, csumvar, csumdist);

			// accumulate our data
			boost::mutex::scoped_lock l(sumLock);
			
			for (int i = 0; i < sumvar.size(); ++i)
			{
				sumvar[i] += csumvar[i];
				sumdist[i] += csumdist[i];
			}
		});
	
	progress p("conservation", (inHits.size() * (inHits.size() + 1)) / 2);
	
	string s(decode(inChain));
	b.put(make_pair(s.c_str(), 0));

	p.step(inHits.size());
	
	for (uint32 i = 0; i + 1 < inHits.size(); ++i)
	{
		b.put(make_pair(inHits[i]->m_aligned.c_str(), i + 1));
		p.step(inHits.size() - i);
	}
	
	b.put(kSentinel);
	threads.join_all();

	for (uint32 i = 0; i < inChain.length(); ++i)
	{
		if (sumdist[i] > 0)
			inResidues[i].m_consweight = sumvar[i] / sumdist[i];
		else
			inResidues[i].m_consweight = 1;
	}

	if (VERBOSE)
		cerr << " done" << endl;
}

// --------------------------------------------------------------------

void CreateHSSP(const MProtein& inProtein, const vector<string>& inDatabanks,
	uint32 inMaxhits, uint32 inMinSeqLength, float inThreshold, float inFragmentCutOff, uint32 inThreads, ostream& inOs)
{
	// construct a set of unique sequences, containing only the largest ones in case of overlap
	vector<sequence> seqset;
	vector<uint32> ix;
	vector<const MChain*> chains;
	string used;
	
	foreach (const MChain* chain, inProtein.GetChains())
	{
		string seq;
		chain->GetSequence(seq);
		
		if (seq.length() < inMinSeqLength)
			continue;
		
		chains.push_back(chain);
		seqset.push_back(encode(seq));
		ix.push_back(ix.size());
	}
	
	if (seqset.empty())
		throw mas_exception(boost::format("Not enough sequences in PDB file of length %1%") % inMinSeqLength);

	if (seqset.size() > 1)
		ClusterSequences(seqset, ix);
	
	// only take the unique sequences
	ix.erase(unique(ix.begin(), ix.end()), ix.end());

	vector<MProfile*> profiles;

	uint32 seqlength = 0;
	foreach (uint32 i, ix)
	{
		const MChain& chain(*chains[ix[i]]);
		
		if (not used.empty())
			used += ", ";
		used += chain.GetChainID();
		
		unique_ptr<MProfile> profile(new MProfile(chain, seqset[ix[i]], inThreshold, inFragmentCutOff));
		
		fs::path blastHits(inProtein.GetID() + '-' + chain.GetChainID() + "-hits.fa");
		fs::ifstream file(blastHits);
		if (not file.is_open())
			throw mas_exception(boost::format("Could not open blast hit file %1%") % blastHits);
		
		{
			progress pr1("processing", fs::file_size(blastHits));
			profile->Process(file, pr1);
		}
		
		CalculateConservation(seqset[ix[i]], profile->m_entries, profile->m_residues);

		seqlength += seqset[ix[i]].length();
		profiles.push_back(profile.release());
	}
	
	stringstream desc;
	if (inProtein.GetHeader().length() >= 50)
		desc << "HEADER     " + inProtein.GetHeader().substr(10, 40) << endl;
	if (inProtein.GetCompound().length() > 10)
		desc << "COMPND     " + inProtein.GetCompound().substr(10) << endl;
	if (inProtein.GetSource().length() > 10)
		desc << "SOURCE     " + inProtein.GetSource().substr(10) << endl;
	if (inProtein.GetAuthor().length() > 10)
		desc << "AUTHOR     " + inProtein.GetAuthor().substr(10) << endl;
	
	CreateHSSPOutput(inProtein.GetID(), desc.str(), inThreshold, inFragmentCutOff, seqlength, chains.size(), ix.size(), used,
		profiles.back()->m_entries, profiles.back()->m_residues, inOs);
}

}

// --------------------------------------------------------------------

#if defined(_MSC_VER)
#include <conio.h>
#endif

int main(int argc, char* argv[])
{
#if P_UNIX
	// enable the dumping of cores to enable postmortem debugging
	rlimit l;
	if (getrlimit(RLIMIT_CORE, &l) == 0)
	{
		l.rlim_cur = l.rlim_max;
		if (l.rlim_cur == 0 or setrlimit(RLIMIT_CORE, &l) < 0)
			cerr << "Failed to set rlimit" << endl;
	}
#endif

	try
	{
		po::options_description desc("MKHSSP options");
		desc.add_options()
			("help,h",							 "Display help message")
			("input,i",		po::value<string>(), "Input PDB file (or PDB ID)")
			("output,o",	po::value<string>(), "Output file, use 'stdout' to output to screen")
			("databank,d",	po::value<vector<string>>(),
												 "Databank(s) to use")
			("threads,a",	po::value<uint32>(), "Number of threads (default is maximum)")
			("min-length",	po::value<uint32>(), "Minimal chain length")
			("max-hits,m",	po::value<uint32>(), "Maximum number of hits to include (default = 1500)")
			("threshold",	po::value<float>(),  "Homology threshold adjustment (default = 0.05)")
			("fragment-cutoff",
							po::value<float>(),  "Minimal alignment length as fraction of chain length (default = 0.75)")
			("verbose,v",						 "Verbose output")
			;
	
		po::positional_options_description p;
		p.add("input", 1);
		p.add("output", 2);
	
		po::variables_map vm;
		po::store(po::command_line_parser(argc, argv).options(desc).positional(p).run(), vm);

		fs::path home = get_home();
		if (fs::exists(home / ".hssprc"))
		{
			fs::ifstream rc(home / ".hssprc");
			po::store(po::parse_config_file(rc, desc), vm);
		}

		po::notify(vm);

		if (vm.count("help") or not vm.count("input") or vm.count("databank") == 0)
		{
			cerr << desc << endl;
			exit(1);
		}

		VERBOSE = vm.count("verbose") > 0;
		
		vector<string> databanks(vm["databank"].as<vector<string>>());
			
		uint32 minlength = 25;
		if (vm.count("min-length"))
			minlength= vm["min-length"].as<uint32>();

		uint32 maxhits = 1500;
		if (vm.count("max-hits"))
			maxhits= vm["max-hits"].as<uint32>();

		float threshold = HSSP::kThreshold;
		if (vm.count("threshold"))
			threshold = vm["threshold"].as<float>();

		float fragmentCutOff = HSSP::kFragmentCutOff;
		if (vm.count("fragment-cutoff"))
			fragmentCutOff = vm["fragment-cutoff"].as<float>();

		uint32 threads = boost::thread::hardware_concurrency();
		if (vm.count("threads"))
			threads = vm["threads"].as<uint32>();
		if (threads < 1)
			threads = 1;
			
		// what input to use
		string input = vm["input"].as<string>();
		io::filtering_stream<io::input> in;
		ifstream infile(input.c_str(), ios_base::in | ios_base::binary);
		if (not infile.is_open())
			throw runtime_error("Error opening input file");

		if (ba::ends_with(input, ".bz2"))
		{
			in.push(io::bzip2_decompressor());
			input.erase(input.length() - 4, string::npos);
		}
		else if (ba::ends_with(input, ".gz"))
		{
			in.push(io::gzip_decompressor());
			input.erase(input.length() - 3, string::npos);
		}
		in.push(infile);

		// Where to write our HSSP file to:
		// either to cout or an (optionally compressed) file.
		ofstream outfile;
		io::filtering_stream<io::output> out;

		if (vm.count("output") and vm["output"].as<string>() != "stdout")
		{
			string output = vm["output"].as<string>();
			outfile.open(output.c_str(), ios_base::out|ios_base::trunc|ios_base::binary);
			
			if (not outfile.is_open())
				throw runtime_error("could not create output file");
			
			if (ba::ends_with(output, ".bz2"))
				out.push(io::bzip2_compressor());
			else if (ba::ends_with(output, ".gz"))
				out.push(io::gzip_compressor());
			out.push(outfile);
		}
		else
			out.push(cout);

		// read protein and calculate the secondary structure
		MProtein a(in);
		a.CalculateSecondaryStructure();
		
		// create the HSSP file
		HSSP::CreateHSSP(a, databanks, maxhits, minlength, threshold, fragmentCutOff, threads, out);
	}
	catch (exception& e)
	{
		cerr << e.what() << endl;
		exit(1);
	}

#if defined(_MSC_VER) && ! NDEBUG
	cerr << "Press any key to quit application ";
	char ch = _getch();
	cerr << endl;
#endif
	
	return 0;
}

