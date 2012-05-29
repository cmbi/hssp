#include "mas.h"

#include <iostream>
#include <set>
#include <cmath>

#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/range/adaptor/sliced.hpp>
#include <boost/filesystem/operations.hpp>
#include <boost/filesystem/convenience.hpp>
#include <boost/filesystem/fstream.hpp>
#include <boost/format.hpp>
#include <boost/foreach.hpp>
#define foreach BOOST_FOREACH
#include <boost/date_time/gregorian/gregorian.hpp>
#include <boost/date_time/date_clock_device.hpp>
#include <boost/regex.hpp>
#include <boost/tr1/cmath.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/device/back_inserter.hpp>
#include <boost/iostreams/copy.hpp>

#include "utils.h"
#include "structure.h"
#include "dssp.h"
#include "matrix.h"
#include "buffer.h"
#include "blast.h"
#include "fetchdbrefs.h"
#include "hssp-nt.h"

using namespace std;
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
	float			m_freq[20];
	float			m_entropy;

	void			Add(uint8 r, float inDistance);
	static MResInfo	NewGap(size_t inDel, float inSumDistance, uint8 inResidue, float inDistance);
	void			AddGap(float inDistance);
};

typedef vector<MResInfo> MResInfoList;

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

MResInfo MResInfo::NewGap(size_t inDel, float inSumDistance, uint8 inResidue, float inDistance)
{
	MResInfo r = {};
	
	r.m_dist[22] = static_cast<uint32>(inDel);
	r.m_sum_dist_weight = r.m_dist_weight[22] = inSumDistance;
	r.Add(inResidue, inDistance);
	
	return r;
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
	
	static MHit*		Create(const string& id, const string& def, const string& seq);

	void				CalculateDistance(const sequence& chain);
	
	string				m_id, m_acc, m_def, m_stid;
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

MHit* MHit::Create(const string& id, const string& def, const string& seq)
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

	void			Process(istream& inHits, float inGapOpen, float inGapExtend, uint32 inMaxHits, uint32 inThreads);
	void			Align(MHit* e, float inGapOpen, float inGapExtend);

	void			AdjustXGapCosts(vector<float>& gop, vector<float>& gep);
	void			AdjustYGapCosts(const sequence& s, vector<float>& gop, vector<float>& gep);

	void			dump(const matrix<float>& B, const matrix<float>& Ix, const matrix<float>& Iy,
						const matrix<int8>& tb, const vector<float>& gopX, const vector<float>& gopY,
						const vector<float>& gepX, const vector<float>& gepY,
						const sequence& sx, const sequence& sy);

	void			PrintStockholm(ostream& os, const string& inChainID, bool inFetchDBRefs) const;
	void			PrintStockholm(ostream& os, const MProtein& inProtein, bool inFetchDBRefs,
						const vector<string>& inUsed, const vector<string>& inAKA) const;

	void			CalculateConservation(uint32 inThreads);

	const MChain&	m_chain;
	sequence		m_seq;
	MResInfoList	m_residues;
	vector<MHit*>	m_entries;
	float			m_threshold, m_frag_cutoff;
	float			m_sum_dist_weight;
};

MProfile::MProfile(const MChain& inChain, const sequence& inSequence, float inThreshold, float inFragmentCutOff)
	: m_chain(inChain), m_seq(inSequence)
	, m_threshold(inThreshold), m_frag_cutoff(inFragmentCutOff), m_sum_dist_weight(0)
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

MProfile::~MProfile()
{
	foreach (auto e, m_entries)
		delete e;
}

void MProfile::dump(const matrix<float>& B, const matrix<float>& Ix, const matrix<float>& Iy,
	const matrix<int8>& tb, const vector<float>& gopX, const vector<float>& gopY,
	const vector<float>& gepX, const vector<float>& gepY, const sequence& sx, const sequence& sy)
{
	ofstream os("alignment.log");
	os.imbue(locale(""));
	
	assert(sx.length() == m_residues.size());
	assert(sx.length() == gopX.size());
	assert(sx.length() == B.dim_m());
	assert(gopY.size() == B.dim_n());
	
	for (uint32 x = 0; x < sx.length(); ++x)
		os << '\t' << (is_gap(sx[x]) ? '.' : kResidues[sx[x]]) << '\t' << gopX[x];
	os << endl;
	
	for (uint32 y = 0; y < sy.length(); ++y)
	{
		os << kResidues[sy[y]];
		for (uint32 x = 0; x < m_residues.size(); ++x)
			os << '\t' << B(x, y) << '\t' << Iy(x, y);
		os << endl
		   << gopY[y];

		for (uint32 x = 0; x < m_residues.size(); ++x)
		{
			switch (tb(x, y))
			{
				case -1:	os << '\t' << Ix(x, y) << '\t' << "|"; break;
				case  0:	os << '\t' << Ix(x, y) << '\t' << "\\"; break;
				case  1:	os << '\t' << Ix(x, y) << '\t' << "-"; break;
				case  2:	os << '\t' << Ix(x, y) << '\t' << "."; break;
			}
		}
		os << endl;
	}
}

void MProfile::AdjustXGapCosts(vector<float>& gop, vector<float>& gep)
{
	assert(gop.size() == m_seq.length());
	assert(gop.size() == m_residues.size());
	
	for (int32 ix = 0; ix < m_residues.size(); ++ix)
	{
		MResInfo& e = m_residues[ix];

		// adjust for secondary structure
		switch (e.m_ss)
		{
			case alphahelix:	gop[ix] *= 2; break;
			case betabridge:	gop[ix] *= 2; break;
			case strand:		break;
			case helix_3:		gop[ix] *= 3; break;
			case helix_5:		gop[ix] *= 2; break;
			case turn:			gop[ix] *= 1.5; break;
			case bend:			break;
			case loop:			break;
		}
		
		// if there is a gap in the alignments, lower gap penalties
		if (e.m_ins > 0)		// gap open penalty is zero when another gap was already created here
			gop[ix] = 0;

		if (e.m_dist[22] > 0)	// lower gap extension penalty for existing gaps here
			gep[ix] = float(e.m_dist[22]) / (m_entries.size() + 1);
		
		// if there is a gap within 8 residues, increase gap penalty
		for (int32 d = 0; d < 8; ++d)
		{
			if (ix + d >= int32(m_residues.size()) or
				m_residues[ix + d].m_dist[22] > 0 or
				ix - d < 0 or
				m_residues[ix - d].m_dist[22] > 0)
			{
				gop[ix] *= (2 + ((8 - d) * 2)) / 8.f;
				break;
			}
		}
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

void MProfile::AdjustYGapCosts(const sequence& s, vector<float>& gop, vector<float>& gep)
{
	for (uint32 y = 0; y < s.length(); ++y)
	{
		if (s[y] < 22)
			gop[y] *= kResidueSpecificPenalty[s[y]];
	}
}

void MProfile::Align(MHit* e, float inGapOpen, float inGapExtend)
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
	
	float logmin = 1.0f / log10(minLength);
	float logdiff = 1.0f + 0.5f * log10(minLength / maxLength);
	
	// initial gap open penalty
	float gop = (inGapOpen / (logdiff * logmin)) * abs(kMPam250MisMatchAverage) * kMPam250ScalingFactor;
	float gep = inGapExtend;

	// position specific gap penalties
	// initial gap extend penalty is adjusted for difference in sequence lengths
	vector<float> gop_a(dimX, gop), gep_a(dimX, gep * (1 + log10(float(dimX) / dimY)));
	AdjustXGapCosts(gop_a, gep_a);
	
	vector<float> gop_b(dimY, gop), gep_b(dimY, gep * (1 + log10(float(dimY) / dimX)));
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
				
				if (highS < M)
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
				Iy(x, y) = max(M - (y < dimY - 1 ? gop_b[y] : 0), Iy1 - gep_b[y]);
			}
			else
			{
				tb(x, y) = -1;
				B(x, y) = Iy1;

				Ix(x, y) = max(M - (x < dimX - 1 ? gop_a[x] : 0), Ix1 - gep_a[x]);
				Iy(x, y) = Iy1 - gep_b[y];
			}
		}
	}

//#if not NDEBUG
//
//bool dmp = false;
//if (e->m_acc == "D2CRJ5")
//{
//	dmp = true;
//	dump(B, Ix, Iy, tb, gop_a, gop_b, gep_a, gep_b, m_seq, e->m_seq);
//}
//
//#endif

	// build the alignment
	x = highX;
	y = highY;

	uint32 ident = 0, similar = 0, length = 0, xgaps = 0;

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
					if (m_seq[x] == '.')
						--length;
					else
						++ident, ++similar;
				}
				else if (score(kMPam250, m_seq[x], e->m_seq[y]) > 0)
					++similar;
				--x;
				--y;
				break;
			
			default:
				assert(false);
				break;
		}
	}
	
	uint32 tix = max(10U, min(length, 80U)) - 10;
	
	// Add the hit only if it is within the required parameters.
	if (length >= m_seq.length() * m_frag_cutoff and				// accept only alignment long enough (suppress fragments)
		ident >= length * (kHomologyThreshold[tix] + m_threshold))	// and those that score high enough
	{
		// reserve space, if needed
		if (xgaps > 0)
		{
			const uint32 kBlockSize = 1024;
			uint32 n = static_cast<uint32>((((m_residues.size() + xgaps) / kBlockSize) + 1) * kBlockSize);
			m_seq.reserve(n);
			foreach (MHit* e, m_entries)
				e->m_aligned.reserve(n);
		}
		
		// update insert/delete counters for the residues
		x = highX;	e->m_ilas = x + 1;
		y = highY;	e->m_jlas = y + 1;

		// trace back to fill aligned sequence and to create gaps in MSA
		e->m_aligned = string(m_seq.length() + xgaps, '.');
		bool gappedx = false, gappedy = false;
		
		while (x >= 0 and y >= 0 and B(x, y) > 0)
		{
			switch (tb(x, y))
			{
				case -1:
					e->m_aligned[x + xgaps] = kResidues[e->m_seq[y]];

					m_residues.insert(m_residues.begin() + x + 1,
						MResInfo::NewGap(m_entries.size() + 1, m_sum_dist_weight, e->m_seq[y], e->m_distance));

					m_seq.insert(m_seq.begin() + x + 1, '.');
					
					foreach (MHit* e, m_entries)
						e->m_aligned.insert(e->m_aligned.begin() + x + 1, '.');
					
					--y;
					--xgaps;
					gappedx = true;
					break;
	
				case 1:
					if (m_residues[x].m_chain_id != 0)
					{
						if (not gappedy)
							++e->m_gaps;
						++e->m_gapn;
						gappedy = true;
					}

					m_residues[x].AddGap(e->m_distance);
					--x;
					break;
	
				case 0:
					e->m_aligned[x + xgaps] = kResidues[e->m_seq[y]];
					m_residues[x].Add(e->m_seq[y], e->m_distance);
					
					if (gappedx or (x + 1 < m_seq.size() and is_gap(m_seq[x + 1]) and not is_gap(m_seq[x])))
						++m_residues[x].m_ins;
					
					if (gappedy)
						++m_residues[x].m_del;
					
					gappedx = gappedy = false;
					
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

		e->m_ifir = x + 2;
		e->m_jfir = y + 2;
		e->m_stid = e->m_acc + '/' +
			boost::lexical_cast<string>(e->m_jfir) + '-' + boost::lexical_cast<string>(e->m_jlas);

		m_entries.push_back(e);

		m_sum_dist_weight += e->m_distance;

//#if not defined(NDEBUG)
//		if (dmp)
//		{
//			string s = decode(m_seq);
//			for (string::size_type i = 72; i < s.length(); i += 73)
//				s.insert(s.begin() + i, '\n');
//			
//			cout << '>' << "PDB" << endl
//				 << s << endl;
//			
//			foreach (MHit* e, m_entries)
//				cout << *e;
//			exit(1);
//		}
//#endif
	}
	else
		delete e;
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

void MProfile::PrintStockholm(ostream& os, const string& inChainID, bool inFetchDBRefs) const
{
	os << "#=GF ID " << inChainID << endl
	   << "#=GF SQ " << m_entries.size() << endl;

	// ## per residue information
	
	int32 nextNr = m_residues.front().m_seq_nr;
	os << "#=GF CC ## RESIDUE INFORMATION" << endl
	   << "#=GF CC SeqNo RESIDUE AA STRUCTURE BP1 BP2  ACC  NOCC VAR" << endl;
	foreach (auto& ri, m_residues)
	{
		if (ri.m_chain_id == 0)
			continue;

		if (ri.m_seq_nr != nextNr)
			os << boost::format("#=GF RI %5.5d       !  !             0   0    0    0    0") % nextNr << endl;

		uint32 ivar = uint32(100 * (1 - ri.m_consweight));
		os << boost::format("#=GF RI %5.5d %s%5.5d  %2.2d") % ri.m_seq_nr % ri.m_dssp % ri.m_nocc % ivar << endl;

		nextNr = ri.m_seq_nr + 1;
	}

	// ## SEQUENCE PROFILE AND ENTROPY
	os << "#=GF CC ## SEQUENCE PROFILE AND ENTROPY" << endl
	   << "#=GF CC   SeqNo PDBNo   V   L   I   M   F   W   Y   G   A   P   S   T   C   H   R   K   Q   E   N   D  NOCC NDEL NINS ENTROPY RELENT WEIGHT" << endl;
	
	nextNr = m_residues.front().m_seq_nr;
	foreach (auto& ri, m_residues)
	{
		if (ri.m_chain_id == 0)
			continue;

		if (ri.m_seq_nr != nextNr)
			os << boost::format("#=GF PR %5.5d           0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0     0    0    0   0.000      0  1.00")
				% nextNr << endl;

		os << boost::format("#=GF PR %5.5d %5.5d %c") % ri.m_seq_nr % ri.m_pdb_nr % ri.m_chain_id;

		for (uint32 i = 0; i < 20; ++i)
			os << boost::format("%4.4d") % uint32(100.0 * ri.m_freq[i] + 0.5);

		uint32 relent = uint32(100 * ri.m_entropy / log(20.0));
		os << "  " << boost::format("%4.4d %4.4d %4.4d   %5.3f   %4.4d  %4.2f") % ri.m_nocc % ri.m_del % ri.m_ins % ri.m_entropy % relent % ri.m_consweight << endl;
		
		nextNr = ri.m_seq_nr + 1;
	}

	// find the longest ID string length
	string::size_type tl = inChainID.length();
	foreach (const MHit* e, m_entries)
	{
		if (tl < e->m_stid.length())
			tl = e->m_stid.length();
	}
	
	os << "#=GS " << inChainID << string(tl - inChainID.length(), ' ') << " CC The query chain" << endl;

	boost::format fmt("#=GS %s HSSP score=%4.2f/%4.2f aligned=%d-%d/%d-%d length=%d ngaps=%d gaplen=%d seqlen=%d");

	foreach (const MHit* e, m_entries)
	{
		string id = e->m_stid + string(tl - e->m_stid.length(), ' ');
		
		os << "#=GS " << id << " ID " << e->m_id << endl
		   << "#=GS " << id << " DE " << e->m_def << endl
		   << fmt % id % e->m_score % (float(e->m_similar) / e->m_length)
				   % e->m_ifir % e->m_ilas % e->m_jfir % e->m_jlas % e->m_length
				   % e->m_gaps % e->m_gapn % e->m_seq.length() << endl;
	
		if (inFetchDBRefs)
		{
			vector<string> pdb;
			const string kBaseURL = "http://mrs.cmbi.ru.nl/mrsws/search/rest/GetLinked/db/uniprot/linkedDatabank/pdb/id/";
			FetchPDBReferences(kBaseURL + e->m_id, pdb);
			if (not pdb.empty())
				os << "#=GS " << id << " DR PDB " << ba::join(pdb, ", ") << endl;
		}
	}

	if (tl < 17)
		tl = 17;

	string::size_type o = 0;
	while (o < m_seq.length())
	{
		string::size_type n = 72;
		if (o + n > m_seq.length())
			n = m_seq.length() - o;
		
		os << endl
		   << inChainID << string(tl - inChainID.length() + 1, ' ') << decode(m_seq.substr(o, n)) << endl;
		
		string ss(n, '.'), ins(n, ' '), del(n, ' '), ent(n, '-'), var(n, '-');
		for (string::size_type i = o; i < o + n; ++i)
		{
			if (m_residues[i].m_seq_nr == 0)
				continue;
			
			switch (m_residues[i].m_ss)
			{
				case alphahelix:	ss[i - o] = 'H'; break;
				case betabridge:	ss[i - o] = 'B'; break;
				case strand:		ss[i - o] = 'E'; break;
				case helix_3:		ss[i - o] = 'G'; break;
				case helix_5:		ss[i - o] = 'I'; break;
				case turn:			ss[i - o] = 'T'; break;
				case bend:			ss[i - o] = 'S'; break;
				case loop:			ss[i - o] = 'C'; break;
			}
			
			ent[i - o] = map_value_to_char(10 * m_residues[i].m_entropy / log(20.0));
			var[i - o] = map_value_to_char(10 * (1 - m_residues[i].m_consweight));
		}
		
		foreach (const MHit* e, m_entries)
			os << e->m_stid << string(tl - e->m_stid.length() + 1, ' ') << e->m_aligned.substr(o, n) << endl;
		
		os << "#=GC SS          " << string(tl - 17 + 1, ' ') << ss << endl
		   << "#=GC Entropy     " << string(tl - 17 + 1, ' ') << ent << endl
		   << "#=GC Variability " << string(tl - 17 + 1, ' ') << var << endl;
		
		o += n;
	}
	
	os << "//" << endl;
}

void MProfile::PrintStockholm(ostream& os, const MProtein& inProtein, bool inFetchDBRefs,
	const vector<string>& inUsed, const vector<string>& inAKA) const
{
	using namespace boost::gregorian;
	date today = day_clock::local_day();

	// write out the profile in Stockholm 1.0 format
	
	os << "# STOCKHOLM 1.0" << endl
	   << "#=GF CC DATE   " << to_iso_extended_string(today) << endl;
	
	string s = inProtein.GetID();
	if (not s.empty())
		os << "#=GF CC PDBID  " << s << endl;
	
	s = inProtein.GetHeader();
	if (not s.empty())
		os << "#=GF CC HEADER " << s << endl;
	
	s = inProtein.GetCompound();
	if (not s.empty())
		os << "#=GF CC COMPND " << s<< endl;
	
	s = inProtein.GetSource();
	if (not s.empty())
		os << "#=GF CC SOURCE " << s << endl;
	
	s = inProtein.GetAuthor();
	if (not s.empty())
		os << "#=GF CC AUTHOR " << s << endl;

	foreach (auto dbref, inProtein.GetDbRef())
		os << "#=GF CC " << dbref << endl;
	
	string queryID = inProtein.GetID();
	if (inProtein.GetChains().size() > 1)
	{
		auto fmt = [](const vector<string>& a) -> string {
			string result;
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
		
		os << "#=GF CC PDB file contains " << inProtein.GetChains().size() << " chains. Used chain" << ( inUsed.size() > 1 ? "s are " : " is " ) << fmt(inUsed) << endl
		   << "#=GF CC Chain " << m_chain.GetChainID() << " is considered to be the same as " << fmt(inAKA) << endl;
		queryID = inProtein.GetID() + '/' + m_chain.GetChainID();
	}
	
	PrintStockholm(os, queryID, inFetchDBRefs);
}

// --------------------------------------------------------------------

void MProfile::Process(istream& inHits, float inGapOpen, float inGapExtend, uint32 inMaxhits, uint32 inThreads)
{
	vector<MHit*> hits;

	string id, def, seq;
	for (;;)
	{
		string line;
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
	
	// sort them by distance
	sort(hits.begin(), hits.end(), [](const MHit* a, const MHit* b) -> bool {
		return a->m_distance < b->m_distance;
	});
	
	// and then align all the hits
	MProgress p2(hits.size(), "aligning");
	foreach (MHit* e, hits)
	{
		Align(e, inGapOpen, inGapExtend);
		p2.Consumed(1);
	}
	
	// sort by score
	sort(m_entries.begin(), m_entries.end(), [](const MHit* a, const MHit* b) -> bool {
		return a->m_score > b->m_score;
	});
	
	if (m_entries.size() > inMaxhits)
		m_entries.erase(m_entries.begin() + inMaxhits, m_entries.end());
	
	CalculateConservation(inThreads);
}

// --------------------------------------------------------------------

// Find the minimal set of overlapping sequences
// In case of strong similarity (distance <= 0.01) we take the longest chain.
void ClusterSequences(vector<sequence>& s, vector<size_t>& ix)
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
// Calculate the variability of a residue, based on dayhoff similarity
// and weights

pair<const char*,uint32> kSentinel((const char*)nullptr, 0);

void CalculateConservation(buffer<pair<const char*,uint32>>& b,
	const vector<MHit*>& inHits, vector<float>& sumvar, vector<float>& sumdist)
{
	size_t length = sumvar.size();
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

void MProfile::CalculateConservation(uint32 inThreads)
{
	vector<float> sumvar(m_seq.length(), 0), sumdist(m_seq.length(), 0);
	
	// Calculate conservation weights in multiple threads to gain speed.
	buffer<pair<const char*,uint32>> b;
	boost::thread_group threads;
	boost::mutex sumLock;

	for (uint32 t = 0; t < inThreads; ++t)
		threads.create_thread([&]() {
			vector<float> csumvar(sumvar.size(), 0), csumdist(sumdist.size(), 0);
			
			HSSP::CalculateConservation(b, m_entries, csumvar, csumdist);

			// accumulate our data
			boost::mutex::scoped_lock l(sumLock);
			
			for (int i = 0; i < sumvar.size(); ++i)
			{
				sumvar[i] += csumvar[i];
				sumdist[i] += csumdist[i];
			}
		});
	
	MProgress p((m_entries.size() * (m_entries.size() + 1)) / 2, "conservation");
	
	string s(decode(m_seq));
	b.put(make_pair(s.c_str(), 0));

	p.Consumed(m_entries.size());
	
	for (uint32 i = 0; i + 1 < m_entries.size(); ++i)
	{
		b.put(make_pair(m_entries[i]->m_aligned.c_str(), i + 1));
		p.Consumed(m_entries.size() - i);
	}
	
	b.put(kSentinel);
	threads.join_all();

	for (uint32 i = 0; i < m_seq.length(); ++i)
	{
		MResInfo& ri = m_residues[i];
		
		if (sumdist[i] > 0)
			ri.m_consweight = sumvar[i] / sumdist[i];
		else
			ri.m_consweight = 1;

		ri.m_entropy = 0;

		static const int8 kResIx[] = {
			//	V   L   I   M   F   W   Y   G   A   P   S   T   C   H   R   K   Q   E   N   D
			   18, 10,  8, 11,  5, 19, 20,  6,  0, 13, 16, 17,  2,  7, 15,  9, 14,  4, 12,  3
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

void CreateHSSP(const MProtein& inProtein, const vector<fs::path>& inDatabanks,
	uint32 inMaxhits, uint32 inMinSeqLength, float inGapOpen, float inGapExtend,
	float inThreshold, float inFragmentCutOff, uint32 inThreads, bool inFetchDBRefs,
	ostream& inOs)
{
	// construct a set of unique sequences, containing only the largest ones in case of overlap
	vector<sequence> seqset;
	vector<size_t> ix;
	vector<const MChain*> chains;
	vector<vector<string>> aka;
	vector<string> used;
	
	foreach (const MChain* chain, inProtein.GetChains())
	{
		string seq;
		chain->GetSequence(seq);
		
		if (seq.length() < inMinSeqLength)
			continue;
		
		chains.push_back(chain);
		seqset.push_back(encode(seq));
		ix.push_back(ix.size());
		aka.push_back(vector<string>());
	}
	
	if (seqset.empty())
		throw mas_exception(boost::format("Not enough sequences in PDB file of length %1%") % inMinSeqLength);

	if (seqset.size() > 1)
		ClusterSequences(seqset, ix);
	
	// only take the unique sequences
	for (size_t i = 0; i < ix.size(); ++i)
	{
		if (ix[i] != i)
			aka[ix[i]].push_back(string(1, chains[i]->GetChainID()));
	}

	ix.erase(unique(ix.begin(), ix.end()), ix.end());

	foreach (size_t i, ix)
		used.push_back(string(1, chains[i]->GetChainID()));

	bool empty = true;

	foreach (size_t i, ix)
	{
		const MChain& chain(*chains[i]);
		
		// do a blast search for inMaxhits * 4 hits.
		vector<char> blastHits;
		blastHits.reserve(inMaxhits * 4 * (chain.GetResidues().size() + 100));
		
		io::filtering_ostream out(io::back_inserter(blastHits));
		SearchAndWriteResultsAsFastA(out, inDatabanks, decode(seqset[i]),
			"blastp", "BLOSUM62", 3, 10, true, true, -1, -1, 4 * inMaxhits,
			inThreads);
		out.flush();

		if (blastHits.empty())
			continue;
		
		MProfile profile(chain, seqset[i], inThreshold, inFragmentCutOff);
		
		io::filtering_istream in(boost::make_iterator_range(blastHits));
		profile.Process(in, inGapOpen, inGapExtend, inMaxhits, inThreads);
		
		if (profile.m_entries.empty())
			continue;

		empty = false;
		profile.PrintStockholm(inOs, inProtein, inFetchDBRefs, used, aka[i]);
	}
	
	if (empty)
		throw mas_exception("No hits found");
}

// --------------------------------------------------------------------

void CreateHSSP(const string& inProtein, const vector<fs::path>& inDatabanks,
	uint32 inMaxhits, uint32 inMinSeqLength, float inGapOpen, float inGapExtend,
	float inThreshold, float inFragmentCutOff, uint32 inThreads, bool inFetchDBRefs,
	ostream& inOs)
{
	MChain* chain = new MChain('A');
	vector<MResidue*>& residues = chain->GetResidues();
	MResidue* last = nullptr;
	uint32 nr = 1;
	foreach (char r, inProtein)
	{
		residues.push_back(new MResidue(nr, r, last));
		++nr;
		last = residues.back();
	}
	
	MProtein protein("INPUT", chain);
	CreateHSSP(protein, inDatabanks, inMaxhits, inMinSeqLength, inGapOpen, inGapExtend,
		inThreshold, inFragmentCutOff, inThreads, inFetchDBRefs, inOs);
}

}
