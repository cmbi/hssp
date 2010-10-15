// 3d dingen

#include "MRS.h"
#include "CSequence.h"

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

#include "matrix.h"
#include "ioseq.h"
#include "utils.h"

#include "CDatabank.h"
#include "CDatabankTable.h"

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

#include <boost/math/quaternion.hpp>

using namespace std;
using namespace tr1;

namespace fs = boost::filesystem;
namespace po = boost::program_options;
namespace io = boost::iostreams;
namespace ba = boost::algorithm;
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

	point&		operator-=(const point& rhs)
				{
					m_x -= rhs.m_x;
					m_y -= rhs.m_y;
					m_z -= rhs.m_z;
					
					return *this;
				}

	point&		operator /=(double f)
				{
					m_x /= f;
					m_y /= f;
					m_z /= f;

					return *this;
				}

	void		rotate(const quaternion& q)
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

point operator-(const point& pt)
{
	return point(-pt.m_x, -pt.m_y, -pt.m_z);
}

double Distance(const point& a, const point& b)
{
	valarray<double> d(3);
	d[0] = a.m_x - b.m_x;
	d[1] = a.m_y - b.m_y;
	d[2] = a.m_z - b.m_z;
	
	d *= d;
	
	return sqrt(d.sum());
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
	
	return os;
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

	if (length > 0.001)
		q /= length;
	else
		q = quaternion(1, 0, 0, 0);

	return q;
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

point centroid(vector<point>& points)
{
	point result;
	
	foreach (point& pt, points)
		result += pt;
	
	result /= points.size();
	
	return result;
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
		
		sum += d.sum();

//cerr << i << ": "
//	 << a[i] << " - " << b[i] << " => "
//	 << d[0] << ',' << d[1] << ',' << d[2]
//	 << " sum " << sum << endl;
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
	int16				mResSeq;
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
						
	void				Rotate(const quaternion& inRotation)
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

	static char MapThreeLetterCode(const string& inThreeLetterCode);

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

char MResidueInfo::MapThreeLetterCode(const string& inThreeLetterCode)
{
	char result = 'U';

	for (uint32 i = 1; i < kResidueTypeCount; ++i)
	{
		if (inThreeLetterCode == kResidueInfo[i].name)
		{
			result = kResidueInfo[i].code;
			break;
		}
	}
	
	return result;
}

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
						
	void				Rotate(const quaternion& inRotation)
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
						
	void				Rotate(const quaternion& inRotation)
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

	point				GetCAlphaPosition(char inChain, int16 inPDBResSeq) const;
	
	void				GetSequence(char inChain, entry& outEntry) const;

	void				GetSequence(char inChain, sequence& outSequence) const;

	void				Center();
	
	void				Translate(const point& inTranslation)
						{
							for (map<char,MChain>::iterator chain = mChains.begin(); chain != mChains.end(); ++chain)
								chain->second.Translate(inTranslation);
						}
						
	void				Rotate(const quaternion& inRotation)
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

point MProtein::GetCAlphaPosition(char inChain, int16 inPDBResSeq) const
{
	if (inChain == 0)
		inChain = mChains.begin()->first;
	
	map<char,MChain>::const_iterator chain = mChains.find(inChain);
	assert(chain != mChains.end());
	
	point result;
	bool found = false;
	
	foreach (const MResidue& r, chain->second.mResidues)
	{
		if (r.mNumber != inPDBResSeq)
			continue;
		
		foreach (const MAtom& a, r.mAtoms)
		{
			if (a.mType == kCarbon and strcmp(a.mName, " CA ") == 0)
			{
				found = true;
				result = a.mLoc;
				break;
			}
		}
	}
	
	if (not found)
		throw mas_exception(boost::format("residue %1% not found in chain %2%") % inPDBResSeq % inChain);
	
	return result;
}

void MProtein::GetSequence(char inChain, entry& outEntry) const
{
	if (inChain == 0)
		inChain = mChains.begin()->first;
	
	map<char,MChain>::const_iterator chain = mChains.find(inChain);
	assert(chain != mChains.end());
	
	string seq;
	
	foreach (const MResidue& r, chain->second.mResidues)
	{
		foreach (const MAtom& a, r.mAtoms)
		{
			if (a.mType == kCarbon and strcmp(a.mName, " CA ") == 0)
			{
				string residueName(a.mResName);
				ba::trim(residueName);
				
				seq += MResidueInfo::MapThreeLetterCode(residueName);
				outEntry.m_positions.push_back(a.mResSeq);
			}
		}
	}
	
	outEntry.m_seq = encode(seq);
}

void MProtein::GetSequence(char inChain, sequence& outSequence) const
{
	if (inChain == 0)
		inChain = mChains.begin()->first;
	
	map<char,MChain>::const_iterator chain = mChains.find(inChain);
	assert(chain != mChains.end());
	
	string seq;
	
	foreach (const MResidue& r, chain->second.mResidues)
	{
		foreach (const MAtom& a, r.mAtoms)
		{
			if (a.mType == kCarbon and strcmp(a.mName, " CA ") == 0)
			{
				string residueName(a.mResName);
				ba::trim(residueName);
				
				seq += MResidueInfo::MapThreeLetterCode(residueName);
			}
		}
	}
	
	outSequence = encode(seq);
}

void MProtein::WritePDB(ostream& os)
{
	for (map<char,MChain>::iterator c = mChains.begin(); c != mChains.end(); ++c)
		c->second.WritePDB(os);
}

double CalculateRMSD(const MProtein& a, const MProtein& b, char chainA, char chainB)
{
	vector<point> pa;
//	a.GetPoints(pa);
	a.GetCAlphaLocations(chainA, pa);
	
	vector<point> pb;
//	b.GetPoints(pb);
	b.GetCAlphaLocations(chainB, pb);
	
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
	bool model = false;
	
	while (not is.eof())
	{
		string line;
		getline(is, line);
		
		if (ba::starts_with(line, "HEADER"))
		{
			prot.mID = line.substr(62, 4);
			continue;
		}
		
		// brain dead support for only the first model in the file
		if (ba::starts_with(line, "MODEL"))
		{
			model = true;
			continue;
		}
		
		if (ba::starts_with(line, "ENDMDL") and model == true)
			break;
		
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
			atom.mResSeq = boost::lexical_cast<int16>(ba::trim_copy(line.substr(22, 4)));
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
			if (line.length() > 76)
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
	matrix<double> M(3, 3, 0);

	for (uint32 i = 0; i < pa.size(); ++i)
	{
		const point& a = pa[i];
		const point& b = pb[i];
		
		M(0, 0) += a.m_x * b.m_x;	M(0, 1) += a.m_x * b.m_y;	M(0, 2) += a.m_x * b.m_z;
		M(1, 0) += a.m_y * b.m_x;	M(1, 1) += a.m_y * b.m_y;	M(1, 2) += a.m_y * b.m_z;
		M(2, 0) += a.m_z * b.m_x;	M(2, 1) += a.m_z * b.m_y;	M(2, 2) += a.m_z * b.m_z;
	}
	
	// Now calculate N, a symmetric 4x4 matrix
	symmetric_matrix<double> N(4);
	
	N(0, 0) =  M(0, 0) + M(1, 1) + M(2, 2);
	N(0, 1) =  M(1, 2) - M(2, 1);
	N(0, 2) =  M(2, 0) - M(0, 2);
	N(0, 3) =  M(0, 1) - M(1, 0);
	
	N(1, 1) =  M(0, 0) - M(1, 1) - M(2, 2);
	N(1, 2) =  M(0, 1) + M(1, 0);
	N(1, 3) =  M(0, 2) + M(2, 0);
	
	N(2, 2) = -M(0, 0) + M(1, 1) - M(2, 2);
	N(2, 3) =  M(1, 2) + M(2, 1);
	
	N(3, 3) = -M(0, 0) - M(1, 1) + M(2, 2);

	// det(N - λI) = 0
	// find the largest λ (λm)
	//
	// Aλ4 + Bλ3 + Cλ2 + Dλ + E = 0
	// A = 1
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
	
	// calculate t = (N - λI)
	matrix<double> li = identity_matrix<double>(4) * lm;
	matrix<double> t = N - li;
	
	// calculate a matrix of cofactors for t
	matrix<double> cf(4, 4);

	const uint32 ixs[4][3] =
	{
		{ 1, 2, 3 },
		{ 0, 2, 3 },
		{ 0, 1, 3 },
		{ 0, 1, 2 }
	};

	uint32 maxR = 0;
	for (uint32 r = 0; r < 4; ++r)
	{
		const uint32* ir = ixs[r];
		
		for (uint32 c = 0; c < 4; ++c)
		{
			const uint32* ic = ixs[c];

			cf(r, c) =
				t(ir[0], ic[0]) * t(ir[1], ic[1]) * t(ir[2], ic[2]) +
				t(ir[0], ic[1]) * t(ir[1], ic[2]) * t(ir[2], ic[0]) +
				t(ir[0], ic[2]) * t(ir[1], ic[0]) * t(ir[2], ic[1]) -
				t(ir[0], ic[2]) * t(ir[1], ic[1]) * t(ir[2], ic[0]) -
				t(ir[0], ic[1]) * t(ir[1], ic[0]) * t(ir[2], ic[2]) -
				t(ir[0], ic[0]) * t(ir[1], ic[2]) * t(ir[2], ic[1]);
		}
		
		if (r > maxR and cf(r, 0) > cf(maxR, 0))
			maxR = r;
	}
	
	// NOTE the negation of the y here, why? Maybe I swapped r/c above?
	quaternion q(cf(maxR, 0), cf(maxR, 1), -cf(maxR, 2), cf(maxR, 3));
	q = normalize(q);
	
	return q;
}

void align_points_iterative(
	MProtein& a, MProtein& b,
	const sequence& sa, const sequence& sb,
	vector<point>& cAlphaA, vector<point>& cAlphaB,
	point& outTranslationA, point& outTranslationB, quaternion& outRotation)
{
	if (cAlphaA.size() != cAlphaB.size())
		throw logic_error("Protein A and B should have the same number of c-alhpa atoms");
	
	if (cAlphaA.size() < 3)
		throw logic_error("Protein A and B should have at least 3 of c-alpha atoms");

	outTranslationA = centroid(cAlphaA);
	foreach (point& pt, cAlphaA)
		pt -= outTranslationA;
	
	outTranslationB = centroid(cAlphaB);
	foreach (point& pt, cAlphaB)
		pt -= outTranslationB;
	
	outRotation = align_points(cAlphaA, cAlphaB);
	foreach (point& pt, cAlphaB)
		pt.rotate(outRotation);

	double angle; point axis;
	tr1::tie(angle, axis) = quaternion_to_angle_axis(outRotation);
	if (VERBOSE > 1)
		cerr << "before iterating, rotation: " << angle << " degrees rotation around axis " << axis << endl;
}

bool align_proteins2(MProtein& a, MProtein& b, char chainA, char chainB,
	point& outTranslationA, point& outTranslationB, quaternion& outRotation)
{
	sequence sa, sb;
	
	a.GetSequence(chainA, sa);
	b.GetSequence(chainB, sb);

	vector<point> cAlphaA, cAlphaB;
	
	a.GetCAlphaLocations(chainA, cAlphaA);
	b.GetCAlphaLocations(chainB, cAlphaB);

	uint32 N = cAlphaA.size();
	uint32 M = cAlphaB.size();
	
	matrix<int32> B(N, M, 0);
	matrix<int8> traceback(N, M, 2);
	int32 high = 0, highA, highB;
	
	for (uint32 ai = 0; ai < cAlphaA.size(); ++ai)
	{
		for (uint32 bi = 0; bi < cAlphaB.size(); ++bi)
		{
			float d = Distance(cAlphaA[ai], cAlphaB[bi]);

			int32 v = 0;
			if (d < 3.5)
				v = 1;

			if (ai > 0 and bi > 0)
				v += B(ai - 1, bi - 1);
			
			int32 ga = -1;
			if (ai > 0)
				ga = B(ai - 1, bi) - 1;
			
			int32 gb = -1;
			if (bi > 0)
				gb = B(ai, bi - 1) - 1;
			
			if (VERBOSE > 5)
				cerr << "d(" << ai << ',' << bi << ") = " << d << " ga: " << ga << " gb: " << gb << " v: " << v << endl;

			if (v >= ga and v >= gb)
			{
				B(ai, bi) = v;
				traceback(ai, bi) = 0;
			}
			else if (ga > gb)
			{
				B(ai, bi) = ga;
				traceback(ai, bi) = 1;
			}
			else
			{
				B(ai, bi) = gb;
				traceback(ai, bi) = -1;
			}
			
			if (B(ai, bi) > high)
			{
				high = B(ai, bi);
				highA = ai;
				highB = bi;
			}
		}
	}

	if (VERBOSE > 4)
		print_matrix(cerr, traceback, sa, sb);
	
	vector<point> newCAlphaA, newCAlphaB;
	
	int32 ai = highA, bi = highB;
	
	while (ai > 0 and bi > 0)
	{
		switch (traceback(ai, bi))
		{
			case 0:
				if (VERBOSE > 3)
					cerr << "map " << ai << '(' << kAA[sa[ai]]
						 << ") to " << bi << '(' << kAA[sb[bi]] << ')' << endl;
				newCAlphaA.push_back(cAlphaA[ai]);
				newCAlphaB.push_back(cAlphaB[bi]);
				--ai;
				--bi;
				break;
			
			case 1:
				--ai;
				break;
			
			case -1:
				--bi;
				break;
			
			default:
				assert(false);
				break;
		}
	}

	outTranslationA = centroid(newCAlphaA);
	foreach (point& pt, newCAlphaA)
		pt -= outTranslationA;
	
	outTranslationB = centroid(newCAlphaB);
	foreach (point& pt, newCAlphaB)
		pt -= outTranslationB;
	
	outRotation = align_points(newCAlphaA, newCAlphaB);

	double angle;
	point axis;

	tr1::tie(angle, axis) = quaternion_to_angle_axis(outRotation);
	if (VERBOSE > 1)
		cerr << "  translation: " << outTranslationA << " and " << outTranslationB << endl
			 << "  rotation: " << angle << " degrees around axis " << axis << endl;

	return angle > 0.01;
}

void align_proteins(MProtein& a, char chainA, MProtein& b, char chainB,
	uint32 iterations,
	substitution_matrix_family& mat, float gop, float gep, float magic)
{
	vector<point> cAlphaA, cAlphaB;

	// fetch sequences
	entry ea(1, a.mID), eb(2, b.mID);
	a.GetSequence(chainA, ea);
	b.GetSequence(chainB, eb);
	
	vector<entry*> aa; aa.push_back(&ea);
	vector<entry*> ab; ab.push_back(&eb);

	vector<entry*> ac;
	joined_node n(new leaf_node(ea), new leaf_node(eb), 0.1, 0.1);

	align(&n, aa, ab, ac, mat, gop, gep, magic, true);

	// now based on this alignment, select c-alpha's that we try to align
	assert(ac.size() == 2);
	assert(ac.front()->m_seq.length() == ac.back()->m_seq.length());
	
	const sequence& sa = ac.front()->m_seq;
	const sequence& sb = ac.back()->m_seq;
	
	const vector<int16>& pa = ac.front()->m_positions;
	const vector<int16>& pb = ac.back()->m_positions;

	sequence nsa, nsb;
	
	for (uint32 i = 0; i < sa.length(); ++i)
	{
		if (sa[i] == kSignalGapCode or sb[i] == kSignalGapCode)
			continue;

		nsa.push_back(sa[i]);
		cAlphaA.push_back(a.GetCAlphaPosition(chainA, pa[i]));

		nsb.push_back(sb[i]);
		cAlphaB.push_back(b.GetCAlphaPosition(chainB, pb[i]));
	}
	
	if (cAlphaA.size() != cAlphaB.size())
		throw logic_error("Protein A and B should have the same number of c-alhpa atoms");
	
	if (cAlphaA.size() < 3)
		throw logic_error("Protein A and B should have at least 3 of c-alpha atoms");
	
	point ta, tb;
	quaternion rotation;
	
	// 
	
	if (cAlphaA.size() != cAlphaB.size())
		throw logic_error("Protein A and B should have the same number of c-alhpa atoms");
	
	if (cAlphaA.size() < 3)
		throw logic_error("Protein A and B should have at least 3 of c-alpha atoms");

	//

	ta = centroid(cAlphaA);
	foreach (point& pt, cAlphaA)
		pt -= ta;
	
	tb = centroid(cAlphaB);
	foreach (point& pt, cAlphaB)
		pt -= tb;
	
	rotation = align_points(cAlphaA, cAlphaB);
	foreach (point& pt, cAlphaB)
		pt.rotate(rotation);

	double angle; point axis;
	tr1::tie(angle, axis) = quaternion_to_angle_axis(rotation);
	if (VERBOSE > 1)
		cerr << "before iterating, rotation: " << angle << " degrees rotation around axis " << axis << endl;
	
	uint32 iteration = 0;
	for (;;)
	{
		a.Translate(-ta);
		b.Translate(-tb);
		b.Rotate(rotation);
		
		if (VERBOSE > 1)
			cerr << "RMSd: (" << iteration << ") " << CalculateRMSD(a, b, chainA, chainB) << endl;

		if (iteration++ >= iterations)
			break;
		
		if (not align_proteins2(a, b, chainA, chainB, ta, tb, rotation))
			break;
	}
}

CDatabankPtr LoadDatabank(
	const string&		inDB)
{
	static CDatabankTable sDBTable;
	return sDBTable.Load(inDB);
}

void align_structures(const string& structureA, const string& structureB,
	uint32 iterations,
	substitution_matrix_family& mat, float gop, float gep, float magic)
{
	MProtein a, b;
	char chainA = 0, chainB = 0;
	
	CDatabankPtr pdb = LoadDatabank("pdb");
	
	string pdbid_a = structureA;
	if (pdbid_a.length() == 5)
	{
		chainA = pdbid_a[4];
		pdbid_a.erase(4, 5);
	}
	else if (pdbid_a.length() != 4)
		throw mas_exception("Please specify a PDB ID in 4 letter code");

	string pdbid_b = structureB;
	if (pdbid_b.length() == 5)
	{
		chainB = pdbid_b[4];
		pdbid_b.erase(4, 5);
	}
	else if (pdbid_b.length() != 4)
		throw mas_exception("Please specify a PDB ID in 4 letter code");

	stringstream file_a(pdb->GetDocument(pdbid_a));
	ParsePDB(file_a, a, false);
	
	if (chainA == 0)
		chainA = a.mChains.begin()->first;

	stringstream file_b(pdb->GetDocument(pdbid_b));
	ParsePDB(file_b, b, false);

	if (chainB == 0)
		chainB = b.mChains.begin()->first;

	align_proteins(a, chainA, b, chainB, iterations, mat, gop, gep, magic);

	if (VERBOSE)
		cerr << "RMSD: " << CalculateRMSD(a, b, chainA, chainB) << endl;
	
	MProtein c;

	c.mChains['A'] = a.mChains[chainA];
	c.mChains['A'].SetChainID('A');

	c.mChains['B'] = b.mChains[chainB];
	c.mChains['B'].SetChainID('B');
	
	ofstream file_o(pdbid_a + chainA + '-' + pdbid_b + chainB + ".pdb");
	c.WritePDB(file_o);
}
