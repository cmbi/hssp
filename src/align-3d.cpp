// 3d dingen

#include "mas.h"

#include <iostream>
#include <iomanip>
#include <string>
#include <limits>
#include <cmath>
#include <numeric>
#include <vector>
#include <map>
#include <valarray>

#include <boost/program_options.hpp>
#include <boost/format.hpp>
#include <boost/filesystem/fstream.hpp>
#include <boost/foreach.hpp>
#define foreach BOOST_FOREACH
#include <boost/tr1/tuple.hpp>
#include <boost/iostreams/device/array.hpp>
#include <boost/iostreams/stream.hpp>
#include <boost/bind.hpp>
#include <boost/algorithm/string.hpp>

#include <boost/array.hpp>

#include <boost/numeric/ublas/symmetric.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/math/quaternion.hpp>

using namespace std;
using namespace tr1;

namespace fs = boost::filesystem;
namespace po = boost::program_options;
namespace io = boost::iostreams;
namespace ba = boost::algorithm;
namespace bu = boost::numeric::ublas;
namespace bm = boost::math;

const double kPI = 4 * std::atan(1);

// --------------------------------------------------------------------

typedef bm::quaternion<double> quaternion;

struct point
{
				point(double x, double y, double z) : m_x(x), m_y(y), m_z(z) {}

				point() : m_x(0), m_y(0), m_z(0) {}

				point(const point& rhs)
					: m_x(rhs.m_x), m_y(rhs.m_y), m_z(rhs.m_z) {}

	point&		operator=(const point& rhs)
				{
					m_x = rhs.m_x;
					m_y = rhs.m_y;
					m_z = rhs.m_z;
					
					return *this;
				}

	point&		operator+=(const point& rhs)
				{
					m_x += rhs.m_x;
					m_y += rhs.m_y;
					m_z += rhs.m_z;
					
					return *this;
				}

	void		rotate(quaternion& q)
				{
					quaternion p(0, m_x, m_y, m_z);
					
					p = q * p * conj(q);

					m_x = p.R_component_2();
					m_y = p.R_component_3();
					m_z = p.R_component_4();
				}

	double		m_x, m_y, m_z;
};

point operator-(const point& lhs, const point& rhs)
{
	return point(lhs.m_x - rhs.m_x, lhs.m_y - rhs.m_y, lhs.m_z - rhs.m_z);
}

ostream& operator<<(ostream& os, const point& pt)
{
	os << '(' << pt.m_x << ',' << pt.m_y << ',' << pt.m_z << ')';
	return os;
}

ostream& operator<<(ostream& os, const vector<point>& pts)
{
	uint32 n = pts.size();
	os << '[' << n << ']';
	
	foreach (const point& pt, pts)
	{
		os << pt;
		if (n-- > 1)
			os << ',';
	}
}

quaternion normalize(quaternion q)
{
	valarray<double> t(4);
	
	t[0] = q.R_component_1();
	t[1] = q.R_component_2();
	t[2] = q.R_component_3();
	t[3] = q.R_component_4();
	
	t *= t;
	
	double length = sqrt(t.sum());
	
	return q / length;
}

tr1::tuple<double,point> quaternion_to_angle_axis(quaternion q)
{
	if (q.R_component_1() > 1)
		q = normalize(q);

	// angle:
	double angle = 2 * acos(q.R_component_1());
	angle = angle * 180 / kPI;

	// axis:
	double s = sqrt(1 - q.R_component_1() * q.R_component_1());
	if (s < 0.001)
		s = 1;
	
	point axis(q.R_component_2() / s, q.R_component_3() / s, q.R_component_4() / s);

	return tr1::make_tuple(angle, axis);
}

point center_points(vector<point>& points)
{
	point t;
	
	foreach (point& pt, points)
	{
		t.m_x += pt.m_x;
		t.m_y += pt.m_y;
		t.m_z += pt.m_z;
	}
	
	t.m_x /= points.size();
	t.m_y /= points.size();
	t.m_z /= points.size();
	
	foreach (point& pt, points)
	{
		pt.m_x -= t.m_x;
		pt.m_y -= t.m_y;
		pt.m_z -= t.m_z;
	}
	
	return t;
}

double rmsd(const vector<point>& a, const vector<point>& b)
{
	double sum = 0;
	for (uint32 i = 0; i < a.size(); ++i)
	{
		valarray<double> d(3);
		
		d[0] = b[i].m_x - a[i].m_x;
		d[1] = b[i].m_y - a[i].m_y;
		d[2] = b[i].m_z - a[i].m_z;
		
		d *= d;
		
		sum += sqrt(d.sum());
	}
	
	return sqrt(sum / a.size());
}

// --------------------------------------------------------------------

struct MResidue;

enum MAtomType
{
	kUnknownAtom,
	kHydrogen,
	// ...
	kCarbon,
	kNitrogen,
	kOxygen,
	kFluorine,
	// ...
	kPhosphorus,
	kSulfur,
	kChlorine,
	kPotassium,
	kCalcium,
	kZinc,
	kSelenium,
	
	kAtomTypeCount
};

struct MAtom
{
	uint32				mSerial;
	char				mName[5];
	char				mAltLoc;
	char				mResName[5];
	char				mChainID;
	uint32				mResSeq;
	char				mICode;
	MAtomType			mType;
	point				mLoc;
	double				mOccupancy;
	double				mTempFactor;
	char				mElement[3];
	int					mCharge;

	void				SetChainID(char inID)
						{
							mChainID = inID;
						}

	void				Translate(const point& inTranslation)
						{
							mLoc += inTranslation;
						}
						
	void				Rotate(quaternion& inRotation)
						{
							mLoc.rotate(inRotation);
						}

	void				WritePDB(ostream& os);
};

void MAtom::WritePDB(ostream& os)
{
	//	1 - 6	Record name "ATOM "
	//	7 - 11	Integer serial Atom serial number.
	//	13 - 16	Atom name Atom name.
	//	17		Character altLoc Alternate location indicator.
	//	18 - 20	Residue name resName Residue name.
	//	22		Character chainID Chain identifier.
	//	23 - 26	Integer resSeq Residue sequence number.
	//	27		AChar iCode Code for insertion of residues.
	//	31 - 38	Real(8.3) x Orthogonal coordinates for X in Angstroms.
	//	39 - 46	Real(8.3) y Orthogonal coordinates for Y in Angstroms.
	//	47 - 54	Real(8.3) z Orthogonal coordinates for Z in Angstroms.
	//	55 - 60	Real(6.2) occupancy Occupancy.
	//	61 - 66	Real(6.2) tempFactor Temperature factor.
	//	77 - 78	LString(2) element Element symbol, right-justified.
	//	79 - 80	LString(2) charge Charge on the atom.
	boost::format atom("ATOM  %5.5d %4.4s%c%3.3s %c%4.4d%c   %8.3f%8.3f%8.3f%6.2f%6.2f          %2.2s%2.2s");

	string charge;
	if (mCharge != 0)
	{
		charge += boost::lexical_cast<string>(abs(mCharge));
		if (mCharge > 0)
			charge += '+';
		else
			charge += '-';
	}
	
	os << (atom % mSerial % mName % mAltLoc % mResName % mChainID % mResSeq % mICode %
		   mLoc.m_x % mLoc.m_y % mLoc.m_z % mOccupancy % mTempFactor % mElement % charge) << endl;
}

enum MResidueType
{
	kUnknownResidue,
	
	//
	kAlanine,				// A	ala
	kArginine,				// R	arg
	kAsparagine,			// N	asn
	kAsparticAcid,			// D	asp
	kCysteine,				// C	cys
	kGlutamicAcid,			// E	glu
	kGlutamine,				// Q	gln
	kGlycine,				// G	gly
	kHistidine,				// H	his
	kIsoleucine,			// I	ile
	kLeucine,				// L	leu
	kLysine,				// K	lys
	kMethionine,			// M	met
	kPhenylalanine,			// F	phe
	kProline,				// P	pro
	kSerine,				// S	ser
	kThreonine,				// T	thr
	kTryptophan,			// W	trp
	kTyrosine,				// Y	tyr
	kValine,				// V	val
	
	kResidueTypeCount
};

struct MResidueInfo
{
	MResidueType		type;
	char				code;
	char				name[4];
} kResidueInfo[] = {
	{ kUnknownResidue,	'U', "UNK" },
	{ kAlanine,			'A', "ALA" },
	{ kArginine,		'R', "ARG" },
	{ kAsparagine,		'N', "ASN" },
	{ kAsparticAcid,	'D', "ASP" },
	{ kCysteine,		'C', "CYS" },
	{ kGlutamicAcid,	'E', "GLU" },
	{ kGlutamine,		'Q', "GLN" },
	{ kGlycine,			'G', "GLY" },
	{ kHistidine,		'H', "HIS" },
	{ kIsoleucine,		'I', "ILE" },
	{ kLeucine,			'L', "LEU" },
	{ kLysine,			'K', "LYS" },
	{ kMethionine,		'M', "MET" },
	{ kPhenylalanine,	'F', "PHE" },
	{ kProline,			'P', "PRO" },
	{ kSerine,			'S', "SER" },
	{ kThreonine,		'T', "THR" },
	{ kTryptophan,		'W', "TRP" },
	{ kTyrosine,		'Y', "TYR" },
	{ kValine,			'V', "VAL" }
};

struct MResidue
{
	int32				mNumber;
	MResidueType		mType;
	vector<MAtom>		mAtoms;

	void				SetChainID(char inID)
						{
							for_each(mAtoms.begin(), mAtoms.end(), boost::bind(&MAtom::SetChainID, _1, inID));
						}

	void				Translate(const point& inTranslation)
						{
							for_each(mAtoms.begin(), mAtoms.end(), boost::bind(&MAtom::Translate, _1, inTranslation));
						}
						
	void				Rotate(quaternion& inRotation)
						{
							for_each(mAtoms.begin(), mAtoms.end(), boost::bind(&MAtom::Rotate, _1, inRotation));
						}

	void				WritePDB(ostream& os)
						{
							for_each(mAtoms.begin(), mAtoms.end(), boost::bind(&MAtom::WritePDB, _1, ref(os)));
						}
};

struct MChain
{
	char				mChainID;
	vector<MResidue>	mResidues;

	void				SetChainID(char inID)
						{
							mChainID = inID;
							for_each(mResidues.begin(), mResidues.end(), boost::bind(&MResidue::SetChainID, _1, inID));
						}

	void				Translate(const point& inTranslation)
						{
							for_each(mResidues.begin(), mResidues.end(), boost::bind(&MResidue::Translate, _1, inTranslation));
						}
						
	void				Rotate(quaternion& inRotation)
						{
							for_each(mResidues.begin(), mResidues.end(), boost::bind(&MResidue::Rotate, _1, inRotation));
						}

	void				WritePDB(ostream& os)
						{
							for_each(mResidues.begin(), mResidues.end(), boost::bind(&MResidue::WritePDB, _1, ref(os)));
							
							boost::format ter("TER    %4.4d      %3.3s %c%4.4d%c");
							
							MResidue& last = mResidues.back();
							
							os << (ter % (last.mAtoms.back().mSerial + 1) % kResidueInfo[last.mType].name % mChainID % last.mNumber % ' ') << endl;
						}
};

struct MProtein
{
	string				mID;
	map<char,MChain>	mChains;
	
	void				GetCAlphaLocations(char inChain, vector<point>& outPoints) const;

	void				Center();
	
//	point				Centroid() const;

	void				Translate(const point& inTranslation)
						{
							for (map<char,MChain>::iterator chain = mChains.begin(); chain != mChains.end(); ++chain)
								chain->second.Translate(inTranslation);
						}
						
	void				Rotate(quaternion& inRotation)
						{
							for (map<char,MChain>::iterator chain = mChains.begin(); chain != mChains.end(); ++chain)
								chain->second.Rotate(inRotation);
						}

	void				WritePDB(ostream& os);
	
	void				GetPoints(vector<point>& outPoints) const
						{
							for (map<char,MChain>::const_iterator chain = mChains.begin(); chain != mChains.end(); ++chain)
							{
								foreach (const MResidue& r, chain->second.mResidues)
								{
									foreach (const MAtom& a, r.mAtoms)
										outPoints.push_back(a.mLoc);
								}
							}
						}
};

void MProtein::Center()
{
	vector<point> p;
	GetPoints(p);
	
	point t = center_points(p);
	
	Translate(point(-t.m_x, -t.m_y, -t.m_z));
}

void MProtein::GetCAlphaLocations(char inChain, vector<point>& outPoints) const
{
	if (inChain == 0)
		inChain = mChains.begin()->first;
	
	map<char,MChain>::const_iterator chain = mChains.find(inChain);
	assert(chain != mChains.end());
	
	foreach (const MResidue& r, chain->second.mResidues)
	{
		foreach (const MAtom& a, r.mAtoms)
		{
			if (a.mType == kCarbon and strcmp(a.mName, " CA ") == 0)
				outPoints.push_back(a.mLoc);
		}
	}
}

void MProtein::WritePDB(ostream& os)
{
	for (map<char,MChain>::iterator c = mChains.begin(); c != mChains.end(); ++c)
		c->second.WritePDB(os);
}

double CalculateRMSD(const MProtein& a, const MProtein& b)
{
	vector<point> pa;
	a.GetPoints(pa);
	
	vector<point> pb;
	b.GetPoints(pb);
	
	return rmsd(pa, pb);
}

inline double ParseFloat(const string& s)
{
	return boost::lexical_cast<double>(ba::trim_copy(s));
}

MResidueType MapResidueName(const string& n)
{
	MResidueType result = kUnknownResidue;
	
	for (uint32 i = 0; i < kResidueTypeCount; ++i)
	{
		if (n == kResidueInfo[i].name)
		{
			result = kResidueInfo[i].type;
			break;
		}
	}
	
	return result;
}

void ParsePDB(istream& is, MProtein& prot, bool cAlhpaOnly)
{
	while (not is.eof())
	{
		string line;
		getline(is, line);
		
		if (ba::starts_with(line, "HEADER"))
		{
			prot.mID = line.substr(63, 4);
			continue;
		}
		
		if (ba::starts_with(line, "ATOM  "))
		{
			if (cAlhpaOnly and line.substr(12, 4) != " CA ")
				continue;
			
			MAtom atom = {};
			//	1 - 6	Record name "ATOM "

			atom.mType = kCarbon;

			atom.mSerial = boost::lexical_cast<uint32>(ba::trim_copy(line.substr(6, 5)));
			//	7 - 11	Integer serial Atom serial number.
			line.copy(atom.mName, 4, 12);
			//	13 - 16	Atom name Atom name.
			atom.mAltLoc = line[16];
			//	17		Character altLoc Alternate location indicator.
			line.copy(atom.mResName, 4, 17);
			//	18 - 20	Residue name resName Residue name.
			atom.mChainID = line[21];
			//	22		Character chainID Chain identifier.
			atom.mResSeq = boost::lexical_cast<uint32>(ba::trim_copy(line.substr(22, 4)));
			//	23 - 26	Integer resSeq Residue sequence number.
			atom.mICode = line[26];
			//	27		AChar iCode Code for insertion of residues.

			atom.mLoc.m_x = ParseFloat(line.substr(30, 8));
			//	31 - 38	Real(8.3) x Orthogonal coordinates for X in Angstroms.
			atom.mLoc.m_y = ParseFloat(line.substr(38, 8));
			//	39 - 46	Real(8.3) y Orthogonal coordinates for Y in Angstroms.
			atom.mLoc.m_z = ParseFloat(line.substr(46, 8));
			//	47 - 54	Real(8.3) z Orthogonal coordinates for Z in Angstroms.
			atom.mOccupancy = ParseFloat(line.substr(54, 6));
			//	55 - 60	Real(6.2) occupancy Occupancy.
			atom.mTempFactor = ParseFloat(line.substr(60, 6));
			//	61 - 66	Real(6.2) tempFactor Temperature factor.
			line.copy(atom.mElement, 2, 76);
			//	77 - 78	LString(2) element Element symbol, right-justified.
			atom.mCharge = 0;
			//	79 - 80	LString(2) charge Charge on the atom.
			
			MResidue residue = {};
			residue.mType = MapResidueName(line.substr(17, 3));
			residue.mNumber = boost::lexical_cast<int32>(ba::trim_copy(line.substr(22, 4)));
			residue.mAtoms.push_back(atom);
			
			char chainID = line[21];
			MChain& chain = prot.mChains[chainID];
			chain.mChainID = chainID;
			chain.mResidues.push_back(residue);
		}
	}
}

// The next function returns the largest solution for a quartic equation
// based on Ferrari's algorithm.
// A depressed quartic is of the form:
//
//   x^4 + ax^2 + bx + c = 0
//
// (since I'm too lazy to find out a better way, I've implemented the
//  routine using complex values to avoid nan's as a result of taking
//  sqrt of a negative number)
double largest_depressed_quartic_solution(double a, double b, double c)
{
	complex<double> P = - (a * a) / 12 - c;
	complex<double> Q = - (a * a * a) / 108 + (a * c) / 3 - (b * b) / 8;
	complex<double> R = - Q / 2.0 + sqrt((Q * Q) / 4.0 + (P * P * P) / 27.0);
	
	complex<double> U = pow(R, 1 / 3.0);
	
	complex<double> y;
	if (U == 0.0)
		y = -5.0 * a / 6.0 + U - pow(Q, 1.0 / 3.0);
	else
		y = -5.0 * a / 6.0 + U - P / (3.0 * U);

	complex<double> W = sqrt(a + 2.0 * y);
	
	// And to get the final result:
	// result = (±W + sqrt(-(3 * alpha + 2 * y ± 2 * beta / W))) / 2;
	// We want the largest result, so:

	valarray<double> t(4);

	t[0] = (( W + sqrt(-(3.0 * a + 2.0 * y + 2.0 * b / W))) / 2.0).real();
	t[1] = (( W + sqrt(-(3.0 * a + 2.0 * y - 2.0 * b / W))) / 2.0).real();
	t[2] = ((-W + sqrt(-(3.0 * a + 2.0 * y + 2.0 * b / W))) / 2.0).real();
	t[3] = ((-W + sqrt(-(3.0 * a + 2.0 * y - 2.0 * b / W))) / 2.0).real();

	return t.max();
}

quaternion align_points(const vector<point>& pa, const vector<point>& pb)
{
	// First calculate M, a 3x3 matrix containing the sums of products of the coordinates of A and B
	
	bu::matrix<double> M(3, 3, 0);

	for (uint32 i = 0; i < pa.size(); ++i)
	{
		const point& a = pa[i];
		const point& b = pb[i];
		
		M(0, 0) += a.m_x * b.m_x;	M(0, 1) += a.m_x * b.m_y;	M(0, 2) += a.m_x * b.m_z;
		M(1, 0) += a.m_y * b.m_x;	M(1, 1) += a.m_y * b.m_y;	M(1, 2) += a.m_y * b.m_z;
		M(2, 0) += a.m_z * b.m_x;	M(2, 1) += a.m_z * b.m_y;	M(2, 2) += a.m_z * b.m_z;
	}

	cerr << "M: " << M << endl;
	
	// Now calculate N, a symmetric 4x4 matrix
	bu::symmetric_matrix<double> N(4, 4);
	
	N(0, 0) = M(0, 0) + M(1, 1) + M(2, 2);
	N(0, 1) = M(1, 2) - M(2, 1);
	N(0, 2) = M(2, 0) - M(0, 2);
	N(0, 3) = M(0, 1) - M(1, 0);
	
	N(1, 1) = M(0, 0) - M(1, 1) - M(2, 2);
	N(1, 2) = M(0, 1) + M(1, 0);
	N(1, 3) = M(0, 2) + M(2, 0);
	
	N(2, 2) = -M(0, 0) + M(1, 1) - M(2, 2);
	N(2, 3) = M(1, 2) + M(2, 1);
	
	N(3, 3) = -M(0, 0) - M(1, 1) + M(2, 2);

	cerr << "N: " << N << endl;

	// det(N - λI) = 0
	// find the largest λ (λm)
	//
	// λ4 + c3λ3 + c2λ2 + c1λ + c0 = 0
	// Aλ4 + Bλ3 + Cλ2 + Dλ + E = 0
	// B = 0
	// and so this is a so-called depressed quartic
	// solve it using Ferrari's algorithm
	
	double C = -2 * (
		M(0, 0) * M(0, 0) + M(0, 1) * M(0, 1) + M(0, 2) * M(0, 2) +
		M(1, 0) * M(1, 0) + M(1, 1) * M(1, 1) + M(1, 2) * M(1, 2) +
		M(2, 0) * M(2, 0) + M(2, 1) * M(2, 1) + M(2, 2) * M(2, 2));
	
	double D = 8 * (M(0, 0) * M(1, 2) * M(2, 1) +
					M(1, 1) * M(2, 0) * M(0, 2) +
					M(2, 2) * M(0, 1) * M(1, 0)) -
			   8 * (M(0, 0) * M(1, 1) * M(2, 2) +
					M(1, 2) * M(2, 0) * M(0, 1) +
					M(2, 1) * M(1, 0) * M(0, 2));
	
	double E = 
		(N(0,0) * N(1,1) - N(0,1) * N(0,1)) * (N(2,2) * N(3,3) - N(2,3) * N(2,3)) +
		(N(0,1) * N(0,2) - N(0,0) * N(2,1)) * (N(2,1) * N(3,3) - N(2,3) * N(1,3)) +
		(N(0,0) * N(1,3) - N(0,1) * N(0,3)) * (N(2,1) * N(2,3) - N(2,2) * N(1,3)) +
		(N(0,1) * N(2,1) - N(1,1) * N(0,2)) * (N(0,2) * N(3,3) - N(2,3) * N(0,3)) +
		(N(1,1) * N(0,3) - N(0,1) * N(1,3)) * (N(0,2) * N(2,3) - N(2,2) * N(0,3)) +
		(N(0,2) * N(1,3) - N(2,1) * N(0,3)) * (N(0,2) * N(1,3) - N(2,1) * N(0,3));
	
	// solve quartic
	double lm = largest_depressed_quartic_solution(C, D, E);
	
	bu::matrix<double> li = bu::identity_matrix<double>(4) * lm;
	bu::matrix<double> t = N - li;
	
	cerr << "t: " << t << endl;
	
	// calculate a matrix of cofactors for t
	bu::matrix<double> cf(4, 4);

	const uint32 ixs[4][3] =
	{
		{ 1, 2, 3 },
		{ 0, 2, 3 },
		{ 0, 1, 3 },
		{ 0, 1, 2 }
	};

	for (uint32 x = 0; x < 4; ++x)
	{
		const uint32* ix = ixs[x];
		
		for (uint32 y = 0; y < 4; ++y)
		{
			const uint32* iy = ixs[y];

			cf(x, y) =
				t(ix[0], iy[0]) * t(ix[1], iy[1]) * t(ix[2], iy[2]) +
				t(ix[0], iy[1]) * t(ix[1], iy[2]) * t(ix[2], iy[0]) +
				t(ix[0], iy[2]) * t(ix[1], iy[0]) * t(ix[2], iy[1]) -
				t(ix[0], iy[2]) * t(ix[1], iy[1]) * t(ix[2], iy[0]) -
				t(ix[0], iy[1]) * t(ix[1], iy[0]) * t(ix[2], iy[2]) -
				t(ix[0], iy[0]) * t(ix[1], iy[2]) * t(ix[2], iy[1]);
		}
	}

	cerr << "cf: " << cf << endl;
	
	// take largest row from this matrix as the quaternion
	double sr = abs(cf(0, 0));
	uint32 lr = 0;
	for (uint32 ri = 1; ri < 4; ++ri)
	{
		double s = abs(cf(ri, 0));
		if (sr < s)
		{
			lr = ri;
			sr = s;
		}
	}
	
	// NOTE the negation of the y here, why?
	quaternion q(cf(lr, 0), cf(lr, 1), -cf(lr, 2), cf(lr, 3));
	q = normalize(q);
	
	return q;
}

quaternion align_proteins(MProtein& a, MProtein& b)
{
	a.Center();
	
	vector<point> cAlphaA;
	a.GetCAlphaLocations(0, cAlphaA);

	b.Center();

	vector<point> cAlphaB;
	b.GetCAlphaLocations(0, cAlphaB);

	if (cAlphaA.size() != cAlphaB.size())
		throw logic_error("Protein A and B should have the same number of c-alhpa atoms");
	
	if (cAlphaA.size() < 3)
		throw logic_error("Protein A and B should have at least 3 of c-alpha atoms");
	
	return align_points(cAlphaA, cAlphaB);
}

void test(double angle)
{
	MProtein a, b;
	
	fs::ifstream file("1crn0.pdb");
	ParsePDB(file, a, false);
	b = a;
	
//	double x = -0.0584427, y = -0.997915, z = -0.0273782;
	double x = 0, y = 0.707107, z = 0.707107;

	angle = kPI * angle / 180.0;
	double s = sin(angle / 2);
	quaternion q1(cos(angle/2), x * s, y * s, z * s);
	
	point axis;
	tr1::tie(angle, axis) = quaternion_to_angle_axis(q1);
	cerr << "q1" << q1 << " is a " << angle << " degrees rotation around axis" << axis << endl;

	b.Rotate(q1);
	
	quaternion q2 = align_proteins(a, b);
	
	tr1::tie(angle, axis) = quaternion_to_angle_axis(q2);
	cerr << "q2" << q2 << " is a " << angle << " degrees rotation around axis" << axis << endl;
	
	b.Rotate(q2);
	
	vector<point> caa, cab;
	a.GetCAlphaLocations(0, caa);
	b.GetCAlphaLocations(0, cab);
	
	cerr << "rmsd: " << rmsd(caa, cab) << endl;
}

int main(int argc, char* argv[])
{
	try
	{
		po::options_description desc("align-3d options");
		desc.add_options()
			("help,h",							 			"Display help message")
			("test,t",										"Simply run test routine")
			("test-angle,a", po::value<double>(),			"Rotate 1crn around this test angle")
			("input,i",		po::value<vector<string> >(),	"Input files")
			;
	
		po::positional_options_description p;
		p.add("input", -1);
	
		po::variables_map vm;
		po::store(po::command_line_parser(argc, argv).options(desc).positional(p).run(), vm);
		po::notify(vm);
		
		if (vm.count("test-angle"))
		{
			test(vm["test-angle"].as<double>());
			exit(0);
		}
		
		if (vm.count("help") or vm.count("input") == 0 or vm["input"].as<vector<string> >().size() != 2)
		{
			cerr << desc << endl;
			exit(1);
		}
		
		MProtein a, b;
		
		fs::ifstream file_a(vm["input"].as<vector<string> >()[0]);
		ParsePDB(file_a, a, false);
		
		fs::ifstream file_b(vm["input"].as<vector<string> >()[1]);
		ParsePDB(file_b, b, false);

		quaternion rotation = align_proteins(a, b);

		double angle;
		point axis;
		tr1::tie(angle, axis) = quaternion_to_angle_axis(rotation);
		cerr << "rotation" << rotation << " is a " << angle << " degrees rotation around axis" << axis << endl;
		
		b.Rotate(rotation);

		cerr << "RMSD: " << CalculateRMSD(a, b) << endl;
		
		MProtein c;

		c.mChains['A'] = a.mChains.begin()->second;
		c.mChains['B'] = b.mChains.begin()->second;
		c.mChains['B'].SetChainID('B');
		
		c.WritePDB(cout);
	}
	catch (const exception& e)
	{
		cerr << e.what() << endl;
		exit(1);
	}
	
	return 0;
}