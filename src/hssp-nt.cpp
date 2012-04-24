#include "mas.h"

#include <iostream>
#include <set>

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
#include <boost/ptr_container/ptr_vector.hpp>
#include <boost/date_time/gregorian/gregorian.hpp>
#include <boost/date_time/date_clock_device.hpp>

#include "utils.h"
#include "structure.h"
#include "dssp.h"

using namespace std;
namespace fs = boost::filesystem;
namespace ba = boost::algorithm;
namespace po = boost::program_options;

int VERBOSE = 0;

// --------------------------------------------------------------------

const float kThreshold = 0.05f;

// precalculated threshold table for identity values between 10 and 80
const double kHomologyThreshold[] = {
	0.795468, 0.75398, 0.717997, 0.686414, 0.658413, 0.633373, 0.610811,
	0.590351, 0.571688, 0.554579, 0.53882, 0.524246, 0.510718, 0.498117,
	0.486344, 0.475314, 0.464951, 0.455194, 0.445984, 0.437275, 0.429023,
	0.421189, 0.413741, 0.406647, 0.399882, 0.39342, 0.38724, 0.381323,
	0.375651, 0.370207, 0.364976, 0.359947, 0.355105, 0.35044, 0.345941,
	0.341599, 0.337406, 0.333352, 0.329431, 0.325636, 0.32196, 0.318396,
	0.314941, 0.311587, 0.308331, 0.305168, 0.302093, 0.299103, 0.296194,
	0.293362, 0.290604, 0.287917, 0.285298, 0.282744, 0.280252, 0.277821,
	0.275448, 0.273129, 0.270865, 0.268652, 0.266488, 0.264372, 0.262302,
	0.260277, 0.258294, 0.256353, 0.254452, 0.252589, 0.250764, 0.248975,
	0.247221
};

bool drop(float score, uint32 length, float threshold)
{
	uint32 ix = max(10U, min(length, 80U)) - 10;
	return score < kHomologyThreshold[ix] + threshold;
}

// --------------------------------------------------------------------
// uBlas compatible matrix types
// matrix is m x n, addressing i,j is 0 <= i < m and 0 <= j < n
// element m i,j is mapped to [i * n + j] and thus storage is row major

template<typename T>
class matrix_base
{
  public:

	typedef T value_type;

	virtual				~matrix_base() {}

	virtual uint32		dim_m() const = 0;
	virtual uint32		dim_n() const = 0;

	virtual value_type&	operator()(uint32 i, uint32 j) { throw std::runtime_error("unimplemented method"); }
	virtual value_type	operator()(uint32 i, uint32 j) const = 0;
	
	matrix_base&		operator*=(const value_type& rhs);

	matrix_base&		operator-=(const value_type& rhs);
};

template<typename T>
matrix_base<T>& matrix_base<T>::operator*=(const T& rhs)
{
	for (uint32 i = 0; i < dim_m(); ++i)
	{
		for (uint32 j = 0; j < dim_n(); ++j)
		{
			operator()(i, j) *= rhs;
		}
	}
	
	return *this;
}

template<typename T>
matrix_base<T>& matrix_base<T>::operator-=(const T& rhs)
{
	for (uint32 i = 0; i < dim_m(); ++i)
	{
		for (uint32 j = 0; j < dim_n(); ++j)
		{
			operator()(i, j) -= rhs;
		}
	}
	
	return *this;
}

template<typename T>
std::ostream& operator<<(std::ostream& lhs, const matrix_base<T>& rhs)
{
	lhs << '[' << rhs.dim_m() << ',' << rhs.dim_n() << ']' << '(';
	for (uint32 i = 0; i < rhs.dim_m(); ++i)
	{
		lhs << '(';
		for (uint32 j = 0; j < rhs.dim_n(); ++j)
		{
			if (j > 0)
				lhs << ',';
			lhs << rhs(i,j);
		}
		lhs << ')';
	}
	lhs << ')';
	
	return lhs;
}

template<typename T>
class matrix : public matrix_base<T>
{
  public:
	typedef T value_type;

						template<typename T2>
						matrix(const matrix_base<T2>& m)
							: m_m(m.dim_m())
							, m_n(m.dim_n())
						{
							m_data = new value_type[m_m * m_n];
							for (uint32 i = 0; i < m_m; ++i)
							{
								for (uint32 j = 0; j < m_n; ++j)
									operator()(i, j) = m(i, j);
							}
						}

						matrix(const matrix& m)
							: m_m(m.m_m)
							, m_n(m.m_n)
						{
							m_data = new value_type[m_m * m_n];
							std::copy(m.m_data, m.m_data + (m_m * m_n), m_data);
						}

	matrix&				operator=(const matrix& m)
						{
							value_type t = new value_type[m.m_m * m.m_n];
							std::copy(m.m_data, m.m_data + (m_m * m_n), t);
							
							delete[] m_data;
							m_data = t;
							m_m = m.m_m;
							m_n = m.m_n;
							
							return *this;
						}
	
						matrix(uint32 m, uint32 n, T v = T())
							: m_m(m)
							, m_n(n)
						{
							m_data = new value_type[m_m * m_n];
							std::fill(m_data, m_data + (m_m * m_n), v);
						}
						
	virtual				~matrix()
						{
							delete [] m_data;
						}
	
	virtual uint32		dim_m() const 					{ return m_m; }
	virtual uint32		dim_n() const					{ return m_n; }

	virtual value_type	operator()(uint32 i, uint32 j) const
						{
							assert(i < m_m); assert(j < m_n);
							return m_data[i * m_n + j];
						}
					
	virtual value_type&	operator()(uint32 i, uint32 j)
						{
							assert(i < m_m); assert(j < m_n);
							return m_data[i * m_n + j];
						}

  private:
	value_type*			m_data;
	uint32				m_m, m_n;
};

// --------------------------------------------------------------------

typedef basic_string<uint8>	sequence;

const int8 kM6Pam250[] = {
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

// 22 real letters and 1 dummy
const char kResidues[] = "ABCDEFGHIKLMNPQRSTVWYZX";

inline uint8 ResidueNr(char inAA)
{
	int result = -1;

	const static uint8 kResidueNrTable[] = {
	//	A   B   C   D   E   F   G   H   I       K   L   M   N       P   Q   R   S   T  U=X  V   W   X   Y   Z
		0,  1,  2,  3,  4,  5,  6,  7,  8, 23,  9, 10, 11, 12, 23, 13, 14, 15, 16, 17, 22, 18, 19, 22, 20, 21
	};
	
	inAA |= 040;
	if (inAA >= 'a' and inAA <= 'z')
		result = kResidueNrTable[inAA - 'a'];
	return result;
}

inline int8 score(const int8 inMatrix[], uint8 inAA1, uint8 inAA2)
{
	int8 result;

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

// --------------------------------------------------------------------

float calculateDistance(const sequence& a, const sequence& b)
{
	const float kDistanceGapOpen = 10;
	const float kDistanceGapExtend = 0.2f;

	int32 x = 0, dimX = a.length();
	int32 y = 0, dimY = b.length();

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
			float M = score(kM6Pam250, a[x], b[y]);
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

struct entry
{
	entry(const entry& e);
	
	entry(const string& id, const string& def, const sequence& seq)
		: m_id(id), m_def(def), m_seq(seq), m_distance(0) {}

	struct insertion
	{
		uint32			m_pos;
		sequence		m_seq;
	};
	
	static entry*		create(const string& id, const string& def,
							const string& seq, const sequence& chain);

	string				m_id, m_def;
	sequence			m_seq, m_aligned;
	float				m_distance;
	vector<insertion>	m_insertions;
};

entry* entry::create(const string& id, const string& def, const string& seq, const sequence& chain)
{
	entry* result = new entry(id, def, encode(seq));
	result->m_aligned = sequence(chain.length(), '-');
	result->m_distance = calculateDistance(result->m_seq, chain);
	return result;
}

ostream& operator<<(ostream& os, const entry& e)
{
	os << e.m_id << '\t' << e.m_distance;
	return os;
}

// --------------------------------------------------------------------

struct residue
{
	uint8			m_chain;
	uint32			m_nocc;
	uint32			m_dist[23];
	float			m_score[23];

	void			add(uint8 r, const int8 m[]);
};

void residue::add(uint8 r, const int8 m[])
{
	assert(r < 23);
	m_nocc += 1;
	m_dist[r] += 1;
	
	for (int i = 0; i < 23; ++i)
	{
		m_score[i] = 0;
		
		for (int j = 0; j < 23; ++j)
			m_score[i] += float(score(m, i, j)) * m_dist[j];
		
		m_score[i] /= m_nocc;
	}
}

struct profile
{
					profile(const sequence& chain, const int8 m[]);

	void			align(entry* e);

	void			dump(ostream& os, const matrix<int8>& tb, const sequence& s);

	sequence		m_chain;
	vector<residue>	m_residues;
	vector<entry*>	m_entries;
};

profile::profile(const sequence& chain, const int8 m[])
	: m_chain(chain)
{
	foreach (uint8 r, chain)
	{
		residue res = { r };
		res.add(r, m);
		m_residues.push_back(res);
	}
}

void profile::dump(ostream& os, const matrix<int8>& tb, const sequence& s)
{
	os << ' ';
	foreach (auto& e, m_residues)
		os << kResidues[e.m_chain];
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

void profile::align(entry* e)
{
	const float kSentinelValue = -(numeric_limits<float>::max() / 2);
	
	int32 x = 0, dimX = m_residues.size();
	int32 y = 0, dimY = e->m_seq.length();
	
#ifdef NDEBUG
	matrix<float> B(dimX, dimY);
	matrix<float> Ix(dimX, dimY);
	matrix<float> Iy(dimX, dimY);
	matrix<int8> tb(dimX, dimY);
#else
	matrix<float> B(dimX, dimY, kSentinelValue);
	matrix<float> Ix(dimX, dimY);
	matrix<float> Iy(dimX, dimY);
	matrix<int8> tb(dimX, dimY, 2);
#endif
	
	const int8* m = kM6Pam250;
	
	float minLength = static_cast<float>(dimX), maxLength = static_cast<float>(dimY);
	if (minLength > maxLength)
		swap(minLength, maxLength);
	
	float gop = 10, gep = 0.2f;
	
	float logmin = 1.0f / log10(minLength);
	float logdiff = 1.0f + 0.5f * log10(minLength / maxLength);
	
	// initial gap open cost, 0.05f is the remaining magical number here...
//	gop = (gop / (logdiff * logmin)) * abs(smat.mismatch_average()) * smat.scale_factor() * magic;

	// position specific gap penalties
	// initial gap extend cost is adjusted for difference in sequence lengths
	vector<float> gop_a(dimX, gop), gep_a(dimX, gep * (1 + log10(float(dimX) / dimY)));
//	adjust_gp(gop_a, gep_a, a);
	
	vector<float> gop_b(dimY, gop), gep_b(dimY, gep * (1 + log10(float(dimY) / dimX)));
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

			float s;
			if (M >= Ix1 and M >= Iy1)
			{
				tb(x, y) = 0;
				B(x, y) = s = M;
			}
			else if (Ix1 >= Iy1)
			{
				tb(x, y) = 1;
				B(x, y) = s = Ix1;
			}
			else
			{
				tb(x, y) = -1;
				B(x, y) = s = Iy1;
			}
			
			if (highS < s)
			{
				highS = s;
				highX = x;
				highY = y;
			}
			
			Ix(x, y) = max(M - (x < dimX - 1 ? gop_a[x] : 0), Ix1 - gep_a[x]);
			Iy(x, y) = max(M - (y < dimY - 1 ? gop_b[y] : 0), Iy1 - gep_b[y]);
		}
	}

	// build the alignment
	x = highX;
	y = highY;

//	if (VERBOSE >= 6)
//		dump(cerr, tb, e->m_seq);

	uint32 ident = 0, length = 0;

	// trace back the matrix
	while (x >= 0 and y >= 0 and B(x, y) > 0)
	{
		++length;
		switch (tb(x, y))
		{
			case -1:
				--y;
				break;

			case 1:
				--x;
				break;

			case 0:
				e->m_aligned[x] = e->m_seq[y];
				if (e->m_seq[y] == m_chain[x])
					++ident;
				--x;
				--y;
				break;
			
			default:
				assert(false);
				break;
		}
	}
	
	if (drop(float(ident) / length, length, kThreshold))
		delete e;
	else
	{
		m_entries.push_back(e);
		for (uint32 i = 0; i < m_chain.length(); ++i)
		{
			if (e->m_aligned[i] == '-')
				continue;
			m_residues[i].add(e->m_aligned[i], m);
		}
	}
}

// --------------------------------------------------------------------

void ProcessHits(istream& data, progress& p, profile& prof)
{
	string id, def, seq;
	for (;;)
	{
		string line;
		getline(data, line);
		if (line.empty() and data.eof())
			break;

		p.step(line.length() + 1);
		
		if (ba::starts_with(line, ">"))
		{
			if (not (id.empty() or seq.empty()))
				prof.align(entry::create(id, def, seq, prof.m_chain));
			
			id.clear();
			def.clear();
			seq.clear();
			
			string::size_type s = line.find(' ');
			id = line.substr(1, s);
			if (s != string::npos)
				def = line.substr(s + 1);
		}
		else
			seq += line;
	}

	if (not (id.empty() or seq.empty()))
		prof.align(entry::create(id, def, seq, prof.m_chain));
}

// --------------------------------------------------------------------

//int main()
//{
//	fs::path infile("tests/1gzm-hits.fa");
//	fs::path outfile("tests/1gzm-a.aln");
//	
////	fs::ifstream file("tests/1crn-hits.fa");
//	fs::ifstream file(infile);
//	fs::ofstream out(outfile);
//	
//	if (file.is_open() and out.is_open())
//	{
////		sequence chain(encode("TTCCPSIVARSNFNVCRLPGTPEAICATYTGCIIIPGATCPGDYAN"));
//
//		sequence chain(encode("MNGTEGPNFYVPFSNKTGVVRSPFEAPQYYLAEPWQFSMLAAYMFLLIMLGFPINFLTLYVTVQHKKLRTPLNYILLNLAVADLFMVFGGFTTTLYTSLHGYFVFGPTGCNLEGFFATLGGEIALWSLVVLAIERYVVVCKPMSNFRFGENHAIMGVAFTWVMALACAAPPLVGWSRYIPEGMQCSCGIDYYTPHEETNNESFVIYMFVVHFIIPLIVIFFCYGQLVFTVKEAAAQQQESATTQKAEKEVTRMVIIMVIAFLICWLPYAGVAFYIFTHQGSDFGPIFMTIPAFFAKTSAVYNPVIYIMMNKQFRNCMVTTLCCGKNDD"));
//		profile p(chain, kM6Pam250);
//				
//		{
//			progress pr1("processing", fs::file_size(infile));
//			ProcessHits(file, pr1, p);
//		}
//		
//		foreach (entry* e, p.m_entries)
//		{
//			string seq;
//			foreach (uint8 r, e->m_aligned)
//			{
//				if (seq.length() % 73 == 72)
//					seq += '\n';
//				seq += r < 23 ? kResidues[r] : '-';
//			}
//			
//			out << '>' << e->m_id << ' ' << e->m_def << endl
//				<< seq << endl;
//		}
//		
////		cout << p << endl;
//	}
//	
//	return 0;
//}


void CreateHSSP(const MProtein& inProtein, const vector<string>& inDatabanks,
	uint32 inMaxhits, uint32 inThreshold, uint32 inThreads, ostream& inOs)
{//
//	uint32 seqlength = 0;
//
//	vector<mseq> alignments(inChainAlignments.size());
//	vector<const MChain*> chains;
//	vector<pair<uint32,uint32> > res_ranges;
//
//	res_list res;
//	hit_list hits;
//
//	uint32 kchain = 0;
//	foreach (string ch, inChainAlignments)
//	{
//		if (ch.length() < 3 or ch[1] != '=')
//			throw mas_exception(boost::format("Invalid chain/stockholm pair specified: '%s'") % ch);
//
//		const MChain& chain = inProtein.GetChain(ch[0]);
//		chains.push_back(&chain);
//
//		string seq;
//		chain.GetSequence(seq);
//
//		// strip off trailing X's. They are not very useful
//		while (ba::ends_with(seq, "X"))
//			seq.erase(seq.end() - 1);
//
//		if (VERBOSE > 1)
//			cerr << "Chain " << ch[0] << " => '" << seq << '\'' << endl;
//
//		seqlength += seq.length();
//		
//		// alignments are stored in datadir
//		fs::path afp = ch.substr(2);
//		if (not fs::exists(afp))
//		{
//			fs::path dataDir = "/data/hssp2/sto/";
//			afp = dataDir / (ch.substr(2) + ".aln.bz2");
//		}
//		if (not fs::exists(afp))
//			throw mas_exception("alignment is missing, exiting");
//
//		fs::ifstream af(afp, ios::binary);
//		if (not af.is_open())
//			throw mas_exception(boost::format("Could not open alignment file '%s'") % afp);
//
//		if (VERBOSE)
//			cerr << "Using fasta file '" << afp << '\'' << endl;
//
//		io::filtering_stream<io::input> in;
//		in.push(io::bzip2_decompressor());
//		in.push(af);
//
//		try {
//			ReadFastA(in, alignments[kchain], seq, inMaxHits, inMinLength, inCutOff);
//		}
//		catch (...)
//		{
//			cerr << "exception while reading file " << afp << endl;
//			throw;
//		}
//
//		// Remove all hits that are not above the threshold here
//		mseq& msa = alignments[kchain];
//		msa.erase(remove_if(msa.begin() + 1, msa.end(), boost::bind(&seq::drop, _1, inCutOff)), msa.end());
//
//		++kchain;
//	}
//
//	string usedChains;
//	kchain = 0;
//	foreach (const MChain* chain, chains)
//	{
//		if (not res.empty())
//			res.push_back(res_ptr(new ResidueHInfo(res.size() + 1)));
//		
//		uint32 first = res.size();
//		
//		mseq& msa = alignments[kchain];
//		ChainToHits(inDatabank, msa, *chain, hits, res);
//		
//		res_ranges.push_back(make_pair(first, res.size()));
//
//		if (not usedChains.empty())
//			usedChains += ',';
//		usedChains += chain->GetChainID();
//
//		++kchain;
//	}
//
//	sort(hits.begin(), hits.end(), compare_hit());
//
//	if (inMaxHits > 0 and hits.size() > inMaxHits)
//		hits.erase(hits.begin() + inMaxHits, hits.end());
//
//	if (hits.empty())
//		throw mas_exception("No hits found or remaining");
//	
//	uint32 nr = 1;
//	foreach (hit_ptr h, hits)
//		h->m_nr = nr++;
//
//	for (uint32 c = 0; c < kchain; ++c)
//	{
//		pair<uint32,uint32> range = res_ranges[c];
//		
//		res_range r(res.begin() + range.first, res.begin() + range.second);
//		CalculateConservation(alignments[c], r);
//
//		foreach (res_ptr ri, r)
//			ri->CalculateVariability(hits);
//	}
//	
//	stringstream desc;
//	if (inProtein.GetHeader().length() >= 50)
//		desc << "HEADER     " + inProtein.GetHeader().substr(10, 40) << endl;
//	if (inProtein.GetCompound().length() > 10)
//		desc << "COMPND     " + inProtein.GetCompound().substr(10) << endl;
//	if (inProtein.GetSource().length() > 10)
//		desc << "SOURCE     " + inProtein.GetSource().substr(10) << endl;
//	if (inProtein.GetAuthor().length() > 10)
//		desc << "AUTHOR     " + inProtein.GetAuthor().substr(10) << endl;
//
//	CreateHSSPOutput(inDatabank, inProtein.GetID(), desc.str(), inCutOff, seqlength,
//		inProtein.GetChains().size(), kchain, usedChains, hits, res, outHSSP);
//	
}

// --------------------------------------------------------------------

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
			("max-hits,m",	po::value<uint32>(), "Maximum number of hits to include (default = 1500)")
			("threshold",	po::value<float>(),  "Homology threshold adjustment (default = 0.05)")
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
			
		uint32 maxhits = 1500;
		if (vm.count("max-hits"))
			maxhits= vm["max-hits"].as<uint32>();

		float threshold = 0.05f;
		if (vm.count("threshold"))
			threshold = vm["threshold"].as<float>();

		uint32 threads = boost::thread::hardware_concurrency();
		if (vm.count("threads"))
			threads = vm["threads"].as<uint32>();
		if (threads < 1)
			threads = 1;
			
		// what input to use
		string input = vm["input"].as<string>();
		io::filtering_stream<io::input> in;
		ifstream infile(input.c_str(), ios_base::in | ios_base::binary);
		if (not input.is_open())
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
		CreateHSSP(a, databanks, maxhits, threshold, threads, out);
	}
	catch (exception& e)
	{
		cerr << e.what() << endl;
		exit(1);
	}

#if P_WIN && P_DEBUG
	cerr << "Press any key to quit application ";
	char ch = _getch();
	cerr << endl;
#endif
	
	return 0;
}

