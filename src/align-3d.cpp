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

const double
	kPI = 4 * std::atan(1),
	kSSBridgeDistance = 3.0,
	kMinimalDistance = 0.5,
	kMinimalCADistance = 9.0,
	kMinHBondEnergy = -9.9,
	kMaxHBondEnergy = -0.5,
	kCouplingConstant = -332 * 0.42 * 0.2,
	kMaxPeptideBondLength = 2.5;

// --------------------------------------------------------------------

inline double ParseFloat(const string& s)
{
	return boost::lexical_cast<double>(ba::trim_copy(s));
}

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

// --------------------------------------------------------------------
// some 3d functions
//
// DihedralAngle

double DotProduct(const point& p1, const point& p2)
{
	return p1.m_x * p2.m_x + p1.m_y * p2.m_y + p1.m_z * p2.m_z;
}

point CrossProduct(const point& p1, const point& p2)
{
	return point(p1.m_y * p2.m_z - p2.m_y * p1.m_z,
				 p1.m_z * p2.m_x - p2.m_z * p1.m_x,
				 p1.m_x * p2.m_y - p2.m_x * p1.m_y);
}

double DihedralAngle(const point& p1, const point& p2, const point& p3, const point& p4)
{
	point v12 = p1 - p2;	// vector from p2 to p1
	point v43 = p4 - p3;	// vector from p3 to p4
	
	point z = p2 - p3;		// vector from p3 to p2
	
	point p = CrossProduct(z, v12);
	point x = CrossProduct(z, v43);
	point y = CrossProduct(z, x);
	
	double u = DotProduct(x, x);
	double v = DotProduct(y, y);
	
	double result = 360;
	if (u > 0 and v > 0)
	{
		u = DotProduct(p, x) / sqrt(u);
		v = DotProduct(p, y) / sqrt(v);
		if (u != 0 or v != 0)
			result = atan2(v, u) * 180 / kPI;
	}
	
	return result;
}

double CosinusAngle(const point& p1, const point& p2, const point& p3, const point& p4)
{
	point v12 = p1 - p2;
	point v34 = p3 - p4;
	
	double result = 0;
	
	double x = DotProduct(v12, v12) * DotProduct(v34, v34);
	if (x > 0)
		result = DotProduct(v12, v34) / sqrt(x);
	
	return result;
}

// --------------------------------------------------------------------

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
	}
	
	return sqrt(sum / a.size());
}

// --------------------------------------------------------------------

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

MAtomType MapElement(const string& inElement)
{
	MAtomType result = kUnknownAtom;
	if (inElement == "H")
		result = kHydrogen;
	else if (inElement == "C")
		result = kCarbon;
	else if (inElement == "N")
		result = kNitrogen;
	else if (inElement == "O")
		result = kOxygen;
	else if (inElement == "F")
		result = kFluorine;
	else if (inElement == "P")
		result = kPhosphorus;
	else if (inElement == "S")
		result = kSulfur;
	else if (inElement == "Cl")
		result = kChlorine;
	else if (inElement == "K")
		result = kPotassium;
	else if (inElement == "Ca")
		result = kCalcium;
	else if (inElement == "Zn")
		result = kZinc;
	else if (inElement == "Se")
		result = kSelenium;
	else
		throw mas_exception(boost::format("Unsupported element %s") % inElement);
	return result;
}

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

	string				GetName() const
						{
							return mName;
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

						operator const point&() const
						{
							return mLoc;
						}

						operator point&()
						{
							return mLoc;
						}
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

class MChain;
class MResidue;

struct MResidueID
{
	char			chain;
	uint16			seqNumber;
	
	bool			operator<(const MResidueID& o) const		{ return chain < o.chain or (chain == o.chain and seqNumber < o.seqNumber); }
};

struct HBond
{
	MResidue*		residue;
	double			energy;
};

enum MBridgeType
{
	nobridge, parallel, antiparallel
};

struct MBridge
{
	char			sheet, ladder;
	MBridgeType		type;
	set<uint16>		link;
	uint16			ib, ie, jb, je, from, to;
};

class MResidue
{
  public:
						MResidue(MChain& chain, uint32 inSeqNumber, uint32 inNumber);

	const MAtom&		GetCAlpha() const;
	const MAtom&		GetC() const;
	const MAtom&		GetN() const;
	const MAtom&		GetO() const;
	const MAtom&		GetH() const				{ return mH; }

	double				Phi() const;
	double				Psi() const;
	tr1::tuple<double,char>
						Alpha() const;
	double				Kappa() const;
	double				TCO() const;
	
	const MResidue*		Next() const;
	const MResidue*		Prev() const;

	void				CheckResidue(MResidue* inPrevious);

	void				SetSSBridgeNr(uint8 inBridgeNr);
	uint8				GetSSBridgeNr() const;

	void				AddAtom(MAtom& inAtom);
	
	HBond*				Donor()						{ return mHBondDonor; }
	HBond*				Acceptor()					{ return mHBondAcceptor; }

	bool				ValidDistance(const MResidue& inNext)
						{
							return Distance(GetC(), inNext.GetN()) <= kMaxPeptideBondLength;
						}

	static bool			TestBond(const MResidue* a, const MResidue* b)
						{
							return a->TestBond(b);
						}

	bool				TestBond(const MResidue* other) const
						{
							return
								(mHBondAcceptor[0].residue == other and mHBondAcceptor[0].energy < kMaxHBondEnergy) or
								(mHBondAcceptor[1].residue == other and mHBondAcceptor[1].energy < kMaxHBondEnergy);
						}

	// bridge functions
	MBridgeType			TestBridge(MResidue* inResidue) const;
	void				ExtendLadder(vector<MBridge>& inBridges);
	void				Sheet(vector<MBridge>& inBridges);
	void				Markstrands(vector<MBridge>& inBridges);

	uint16				GetSeqNumber() const
						{
							return mSeqNumber;
						}

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

	void				WriteDSSP(ostream& os);

	static double		CalculateHBondEnergy(MResidue& inDonor, MResidue& inAcceptor);

//  protected:
	MChain&				mChain;
	int32				mSeqNumber, mNumber;
	MResidueType		mType;
	uint8				mSSBridgeNr;
	MAtom				mH;
	HBond				mHBondDonor[2], mHBondAcceptor[2];
	vector<MAtom>		mAtoms;
};

struct MChain
{
	char				mChainID;
	vector<MResidue*>	mResidues;

	void				SetChainID(char inID)
						{
							mChainID = inID;
							for_each(mResidues.begin(), mResidues.end(), boost::bind(&MResidue::SetChainID, _1, inID));
						}

	void				CheckResidues()
						{
							if (not mResidues.empty())
								mResidues.front()->CheckResidue(nil);

							for (uint32 i = 1; i < mResidues.size(); ++i)
								mResidues[i]->CheckResidue(mResidues[i - 1]);
						}
						
	MResidue&			GetResidueBySeqNumber(uint16 inSeqNumber);

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
							
							MResidue& last = *mResidues.back();
							
							os << (ter % (last.mAtoms.back().mSerial + 1) % kResidueInfo[last.mType].name % mChainID % last.mNumber % ' ') << endl;
						}

	void				WriteDSSP(ostream& os)
						{
							for_each(mResidues.begin(), mResidues.end(), boost::bind(&MResidue::WriteDSSP, _1, ref(os)));
						}

	const MResidue*		Next(const MResidue* res) const;
	const MResidue*		Prev(const MResidue* res) const;

};

struct MProtein
{
						MProtein() {}
						
						MProtein(istream& is, bool inCAlphaOnly = false);

	void				CalculateSecondaryStructure();
	
	void				GetCAlphaLocations(char inChain, vector<point>& outPoints) const;

	point				GetCAlphaPosition(char inChain, int16 inPDBResSeq) const;
	
	void				GetSequence(char inChain, entry& outEntry) const;

	void				GetSequence(char inChain, sequence& outSequence) const;
	
	void				CalculateSSBridges();

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

	void				WriteDSSP(ostream& os);
	
	void				GetPoints(vector<point>& outPoints) const
						{
							for (map<char,MChain>::const_iterator chain = mChains.begin(); chain != mChains.end(); ++chain)
							{
								foreach (const MResidue* r, chain->second.mResidues)
								{
									foreach (const MAtom& a, r->mAtoms)
										outPoints.push_back(a.mLoc);
								}
							}
						}

	MResidue&			GetResidue(char inChainID, uint16 inSeqNumber);

//  private:

	string				mID;
	map<char,MChain>	mChains;
	
	vector<pair<MResidueID,MResidueID>>
						mSSBonds;
};

// --------------------------------------------------------------------

MResidue::MResidue(MChain& chain, uint32 inSeqNumber, uint32 inNumber)
	: mChain(chain)
	, mSeqNumber(inSeqNumber)
	, mNumber(inNumber)
	, mType(kUnknownResidue)
	, mSSBridgeNr(0)
{
	mHBondDonor[0].energy = mHBondDonor[1].energy = mHBondAcceptor[0].energy = mHBondAcceptor[1].energy = 0;
	mHBondDonor[0].residue = mHBondDonor[1].residue = mHBondAcceptor[0].residue = mHBondAcceptor[1].residue = nil;
}

void MResidue::CheckResidue(MResidue* inPrevious)
{
	mH = GetN();
	
	if (mType != kProline and inPrevious != nil)
	{
		const MAtom& pc = inPrevious->GetC();
		const MAtom& po = inPrevious->GetO();
		
		double CODistance = Distance(pc, po);
		
		mH.mLoc.m_x += (pc.mLoc.m_x - po.mLoc.m_x) / CODistance; 
		mH.mLoc.m_y += (pc.mLoc.m_y - po.mLoc.m_y) / CODistance; 
		mH.mLoc.m_z += (pc.mLoc.m_z - po.mLoc.m_z) / CODistance; 
	}
}

void MResidue::AddAtom(MAtom& inAtom)
{
	if (inAtom.mType == kHydrogen)
	{
	}
	else
		mAtoms.push_back(inAtom);
}

const MAtom& MResidue::GetCAlpha() const
{
	vector<MAtom>::const_iterator atom = find_if(mAtoms.begin(), mAtoms.end(),
		boost::bind(&MAtom::GetName, _1) == " CA ");
	if (atom == mAtoms.end())
		throw mas_exception(boost::format("Residue %d has no C-alpha atom") % mNumber);
	return *atom;
}

const MAtom& MResidue::GetC() const
{
	vector<MAtom>::const_iterator atom = find_if(mAtoms.begin(), mAtoms.end(),
		boost::bind(&MAtom::GetName, _1) == " C  ");
	if (atom == mAtoms.end())
		throw mas_exception(boost::format("Residue %d has no c atom") % mNumber);
	return *atom;
}

const MAtom& MResidue::GetN() const
{
	vector<MAtom>::const_iterator atom = find_if(mAtoms.begin(), mAtoms.end(),
		boost::bind(&MAtom::GetName, _1) == " N  ");
	if (atom == mAtoms.end())
		throw mas_exception(boost::format("Residue %d has no N atom") % mNumber);
	return *atom;
}

const MAtom& MResidue::GetO() const
{
	vector<MAtom>::const_iterator atom = find_if(mAtoms.begin(), mAtoms.end(),
		boost::bind(&MAtom::GetName, _1) == " O  ");
	if (atom == mAtoms.end())
		throw mas_exception(boost::format("Residue %d has no O atom") % mNumber);
	return *atom;
}

const MResidue* MResidue::Next() const
{
	return mChain.Next(this);
}

const MResidue* MResidue::Prev() const
{
	return mChain.Prev(this);
}

double MResidue::Phi() const
{
	double result = 360;
	const MResidue* prev = Prev();
	if (prev != nil)
		result = DihedralAngle(prev->GetC(), GetN(), GetCAlpha(), GetC());
	return result;
}

double MResidue::Psi() const
{
	double result = 360;
	const MResidue* next = Next();
	if (next != nil)
		result = DihedralAngle(GetN(), GetCAlpha(), GetC(), next->GetN());
	return result;
}

tr1::tuple<double,char> MResidue::Alpha() const
{
	double alhpa = 360;
	char chirality = ' ';
	
	const MResidue* prev = Prev();
	const MResidue* next = Next();
	const MResidue* nextNext = next ? next->Next() : nil;
	if (prev != nil and nextNext != nil)
	{
		alhpa = DihedralAngle(prev->GetCAlpha(), GetCAlpha(), next->GetCAlpha(), nextNext->GetCAlpha());
		if (alhpa < 0)
			chirality = '-';
		else
			chirality = '+';
	}
	return tr1::make_tuple(alhpa, chirality);
}

double MResidue::Kappa() const
{
	double result = 360;
	const MResidue* prev = Prev();
	const MResidue* prevPrev = prev ? prev->Prev() : nil;
	const MResidue* next = Next();
	const MResidue* nextNext = next ? next->Next() : nil;
	if (prevPrev != nil and nextNext != nil)
	{
		double ckap = CosinusAngle(GetCAlpha(), prevPrev->GetCAlpha(), nextNext->GetCAlpha(), GetCAlpha());
		double skap = sqrt(1 - ckap * ckap);
		result = atan2(skap, ckap) * 180 / kPI;
	}
	return result;
}

double MResidue::TCO() const
{
	double result = 0;
	const MResidue* prev = Prev();
	if (prev != nil)
		result = CosinusAngle(GetC(), GetO(), prev->GetC(), prev->GetO());
	return result;
}

void MResidue::SetSSBridgeNr(uint8 inBridgeNr)
{
	if (mType != kCysteine)
		throw mas_exception("Only cysteine residues can form sulphur bridges");
	mSSBridgeNr = inBridgeNr;
}

uint8 MResidue::GetSSBridgeNr() const
{
	if (mType != kCysteine)
		throw mas_exception("Only cysteine residues can form sulphur bridges");
	return mSSBridgeNr;
}

double MResidue::CalculateHBondEnergy(MResidue& inDonor, MResidue& inAcceptor)
{
	double result = 0;
	
	if (inDonor.mType != kProline)
	{
		double distanceHO = Distance(inDonor.GetH(), inAcceptor.GetO());
		double distanceHC = Distance(inDonor.GetH(), inAcceptor.GetC());
		double distanceNC = Distance(inDonor.GetN(), inAcceptor.GetC());
		double distanceNO = Distance(inDonor.GetN(), inAcceptor.GetO());
		
		if (distanceHO < kMinimalDistance or distanceHC < kMinimalDistance or distanceNC < kMinimalDistance or distanceNO < kMinimalDistance)
			result = kMinHBondEnergy;
		else
			result = kCouplingConstant / distanceHO - kCouplingConstant / distanceHC + kCouplingConstant / distanceNC - kCouplingConstant / distanceNO;

		if (result < kMinHBondEnergy)
			result = kMinHBondEnergy;
	}

	// update donor
	if (result < inDonor.mHBondAcceptor[0].energy)
	{
		inDonor.mHBondAcceptor[1] = inDonor.mHBondAcceptor[0];
		inDonor.mHBondAcceptor[0].residue = &inAcceptor;
		inDonor.mHBondAcceptor[0].energy = result;
	}
	else if (result < inDonor.mHBondAcceptor[1].energy)
	{
		inDonor.mHBondAcceptor[1].residue = &inAcceptor;
		inDonor.mHBondAcceptor[1].energy = result;
	}		

	// and acceptor
	if (result < inAcceptor.mHBondDonor[0].energy)
	{
		inAcceptor.mHBondDonor[1] = inAcceptor.mHBondDonor[0];
		inAcceptor.mHBondDonor[0].residue = &inDonor;
		inAcceptor.mHBondDonor[0].energy = result;
	}
	else if (result < inAcceptor.mHBondDonor[1].energy)
	{
		inAcceptor.mHBondDonor[1].residue = &inDonor;
		inAcceptor.mHBondDonor[1].energy = result;
	}		
	
	return result;
}

MBridgeType MResidue::TestBridge(MResidue* test) const
{
	const MResidue* i0 = Prev();
	const MResidue* i1 = this;
	const MResidue* i2 = Next();
	const MResidue* j0 = test->Prev();
	const MResidue* j1 = test;
	const MResidue* j2 = test->Next();
	
	MBridgeType result = nobridge;
	if (i0 and i2 and j0 and j2)
	{
		if ((TestBond(i2, j1) and TestBond(j1, i0)) or (TestBond(j2, i1) and TestBond(i1, j0)))
			result = parallel;
		else if ((TestBond(i2, j0) and TestBond(j2, i0)) or (TestBond(j1, i1) and TestBond(i1, j1)))
			result = antiparallel;
	}
	
	return result;
}

void MResidue::ExtendLadder(vector<MBridge>& inBridges)
{
}

void MResidue::Sheet(vector<MBridge>& inBridges)
{
}

void MResidue::Markstrands(vector<MBridge>& inBridges)
{
}

void MResidue::WriteDSSP(ostream& os)
{
/*   
  #  RESIDUE AA STRUCTURE BP1 BP2  ACC     N-H-->O    O-->H-N    N-H-->O    O-->H-N    TCO  KAPPA ALPHA  PHI   PSI    X-CA   Y-CA   Z-CA 
 */
	boost::format kDSSPResidueLine(
	"%5.5d%5.5d %c %c  %c.....%c.....           %11s%11s%11s%11s  %6.3f%6.1f%6.1f%6.1f%6.1f %6.1f %6.1f %6.1f");
	
	const MAtom& ca = GetCAlpha();
	
	char code = kResidueInfo[mType].code;
	if (mType == kCysteine and GetSSBridgeNr() != 0)
		code = 'a' - 1 + (GetSSBridgeNr() % 26);

	double alpha;
	char chirality;
	tr1::tie(alpha,chirality) = Alpha();
	
	string NHO1 = "0, 0.0", NHO2 = "0, 0.0", ONH1 = "0, 0.0", ONH2 = "0, 0.0";
	
	MResidue* acceptor = mHBondAcceptor[0].residue;
	if (acceptor != nil)
	{
		int32 d = acceptor->mNumber - mNumber;
		NHO1 = (boost::format("%d,%3.1f") % d % mHBondAcceptor[0].energy).str();
	}
	
	acceptor = mHBondAcceptor[1].residue;
	if (acceptor != nil)
	{
		int32 d = acceptor->mNumber - mNumber;
		NHO2 = (boost::format("%d,%3.1f") % d % mHBondAcceptor[1].energy).str();
	}

	MResidue* donor = mHBondDonor[0].residue;
	if (donor != nil)
	{
		int32 d = donor->mNumber - mNumber;
		ONH1 = (boost::format("%d,%3.1f") % d % mHBondDonor[0].energy).str();
	}
	
	donor = mHBondDonor[1].residue;
	if (donor != nil)
	{
		int32 d = donor->mNumber - mNumber;
		ONH2 = (boost::format("%d,%3.1f") % d % mHBondDonor[1].energy).str();
	}
		
	cout << (kDSSPResidueLine % mNumber % ca.mResSeq % ca.mChainID % code % ' ' % chirality % 
		NHO1 % ONH1 % NHO2 % ONH2 %
		TCO() % Kappa() % alpha % Phi() % Psi() % ca.mLoc.m_x % ca.mLoc.m_y % ca.mLoc.m_z) << endl;
}

// --------------------------------------------------------------------

MResidue& MChain::GetResidueBySeqNumber(uint16 inSeqNumber)
{
	vector<MResidue*>::iterator r = find_if(mResidues.begin(), mResidues.end(),
		boost::bind(&MResidue::GetSeqNumber, _1) == inSeqNumber);
	if (r == mResidues.end())
		throw mas_exception(boost::format("Residue %d not found") % inSeqNumber);
	return **r;
}

const MResidue* MChain::Next(const MResidue* res) const
{
	const MResidue* result = nil;
	for (uint32 i = 0; i < mResidues.size() - 1; ++i)
	{
		if (mResidues[i] == res)
		{
			result = mResidues[i + 1];
			break;
		}
	}
	return result;
}

const MResidue* MChain::Prev(const MResidue* res) const
{
	const MResidue* result = nil;
	for (uint32 i = 1; i < mResidues.size(); ++i)
	{
		if (mResidues[i] == res)
		{
			result = mResidues[i - 1];
			break;
		}
	}
	return result;
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

// --------------------------------------------------------------------

MProtein::MProtein(istream& is, bool cAlphaOnly)
{
	bool model = false;
	MResidue* residue = NULL;
	uint32 resNumber = 0;
	
	while (not is.eof())
	{
		string line;
		getline(is, line);
		
		if (ba::starts_with(line, "HEADER"))
		{
			mID = line.substr(62, 4);
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
		
		if (ba::starts_with(line, "SSBOND"))
		{
			//SSBOND   1 CYS A    6    CYS A   11                          1555   1555  2.03  
			pair<MResidueID,MResidueID> ssbond;
			ssbond.first.chain = line[15];
			ssbond.first.seqNumber = boost::lexical_cast<uint16>(ba::trim_copy(line.substr(16, 5)));
			ssbond.second.chain = line[29];
			ssbond.second.seqNumber = boost::lexical_cast<uint16>(ba::trim_copy(line.substr(30, 5)));

			mSSBonds.push_back(ssbond);
			continue;
		}
		
		if (ba::starts_with(line, "ATOM  "))
		{
			if (cAlphaOnly and line.substr(12, 4) != " CA ")
				continue;
			
			MAtom atom = {};
			//	1 - 6	Record name "ATOM "

			//	7 - 11	Integer serial Atom serial number.
			atom.mSerial = boost::lexical_cast<uint32>(ba::trim_copy(line.substr(6, 5)));
			//	13 - 16	Atom name Atom name.
			line.copy(atom.mName, 4, 12);
			//	17		Character altLoc Alternate location indicator.
			atom.mAltLoc = line[16];
			//	18 - 20	Residue name resName Residue name.
			line.copy(atom.mResName, 4, 17);
			//	22		Character chainID Chain identifier.
			atom.mChainID = line[21];
			//	23 - 26	Integer resSeq Residue sequence number.
			atom.mResSeq = boost::lexical_cast<int16>(ba::trim_copy(line.substr(22, 4)));
			//	27		AChar iCode Code for insertion of residues.
			atom.mICode = line[26];

			//	31 - 38	Real(8.3) x Orthogonal coordinates for X in Angstroms.
			atom.mLoc.m_x = ParseFloat(line.substr(30, 8));
			//	39 - 46	Real(8.3) y Orthogonal coordinates for Y in Angstroms.
			atom.mLoc.m_y = ParseFloat(line.substr(38, 8));
			//	47 - 54	Real(8.3) z Orthogonal coordinates for Z in Angstroms.
			atom.mLoc.m_z = ParseFloat(line.substr(46, 8));
			//	55 - 60	Real(6.2) occupancy Occupancy.
			atom.mOccupancy = ParseFloat(line.substr(54, 6));
			//	61 - 66	Real(6.2) tempFactor Temperature factor.
			atom.mTempFactor = ParseFloat(line.substr(60, 6));
			//	77 - 78	LString(2) element Element symbol, right-justified.
			if (line.length() > 76)
				line.copy(atom.mElement, 2, 76);
			//	79 - 80	LString(2) charge Charge on the atom.
			atom.mCharge = 0;
			
			atom.mType = MapElement(ba::trim_copy(line.substr(77, 2)));
			
			char chainID = line[21];
			MChain& chain = mChains[chainID];
			chain.mChainID = chainID;

			if (residue == NULL or residue->mSeqNumber != atom.mResSeq)
			{
				if (chain.mResidues.empty())
					++resNumber;
				
				residue = new MResidue(chain, atom.mResSeq, resNumber);
				++resNumber;

				residue->mType = MapResidueName(line.substr(17, 3));

				chain.mResidues.push_back(residue);
			}
				
			residue->AddAtom(atom);
		}
	}
}

void MProtein::CalculateSecondaryStructure()
{
	// validate all residues and calculate the location of the backbone Hydrogen
	for (map<char,MChain>::iterator chain = mChains.begin(); chain != mChains.end(); ++chain)
		chain->second.CheckResidues();
	
	// map the sulfur bridges
	sort(mSSBonds.begin(), mSSBonds.end());
	
	typedef pair<MResidueID,MResidueID> SSBond;
	
	uint32 ssbondNr = 1;
	foreach (const SSBond& ssbond, mSSBonds)
	{
		MResidue& first = GetResidue(ssbond.first.chain, ssbond.first.seqNumber);
		MResidue& second = GetResidue(ssbond.second.chain, ssbond.second.seqNumber);
		
		first.SetSSBridgeNr(ssbondNr);
		second.SetSSBridgeNr(ssbondNr);
		++ssbondNr;
	}
	
	// Calculate the HBond energies
	vector<MResidueID> residues;
	for (map<char,MChain>::iterator chain = mChains.begin(); chain != mChains.end(); ++chain)
	{
		MResidueID id = { chain->first };
		foreach (const MResidue* r, chain->second.mResidues)
		{
			id.seqNumber = r->mSeqNumber;
			residues.push_back(id);
		}
	}
	
	for (uint32 i = 0; i < residues.size() - 1; ++i)
	{
		MResidue& ri = GetResidue(residues[i].chain, residues[i].seqNumber);
		
		for (uint32 j = i + 1; j < residues.size(); ++j)
		{
			MResidue& rj = GetResidue(residues[j].chain, residues[j].seqNumber);
			
			if (Distance(ri.GetCAlpha(), rj.GetCAlpha()) < kMinimalCADistance)
			{
				MResidue::CalculateHBondEnergy(ri, rj);
				if (j != i + 1)
					MResidue::CalculateHBondEnergy(rj, ri);
			}
		}
	}
	
	// Calculate Bridges
	vector<MBridge> bridges;
	for (uint32 i = 1; i < residues.size() - 4; ++i)
	{
		MResidue& ir = GetResidue(residues[i].chain, residues[i].seqNumber);
		
		for (uint32 j = i + 3; j < residues.size() - 1; ++j)
		{
			MResidue& jr = GetResidue(residues[j].chain, residues[j].seqNumber);
			
			MBridgeType type = ir.TestBridge(&jr);
			if (type != nobridge)
			{
				bool found = false;
				foreach (MBridge& bridge, bridges)
				{
					if (type == parallel and bridge.type == parallel and
						bridge.ie + 1 == i and bridge.je + 1 == j and
						residues[i - 1].chain == residues[i].chain and
						residues[j - 1].chain == residues[j].chain)
					{
						++bridge.ie;
						++bridge.je;
						found = true;
						break;
					}

					if (type == antiparallel and bridge.type == antiparallel and
						bridge.ie + 1 == i and bridge.jb - 1 == j and
						residues[i - 1].chain == residues[i].chain and
						residues[j + 1].chain == residues[j].chain)
					{
						++bridge.ie;
						--bridge.jb;
						found = true;
						break;
					}
				}
				
				if (not found)
				{
					MBridge bridge = {};
					
					bridge.ib = bridge.ie = i;
					bridge.jb = bridge.je = j;
					bridges.push_back(bridge);
				}
			}
		}
	}

	if (not bridges.empty())
	{
		// extend ladders
		for (uint32 i = 0; i < bridges.size(); ++i)
		{
			MBridge& bridge = bridges[i];
			
			for (uint32 j = i + 1; j < bridges.size() and bridge.to == 0; ++j)
			{
				uint32 ib1 = bridges[j].ib;
				uint32 jb1 = bridges[j].jb;
				uint32 je1 = bridges[j].je;

				bool bulge = residues[bridge.ie].chain == residues[ib1].chain and
					ib1 - bridge.ie < 6 and
					bridges[j].type == bridge.type and
					bridges[j].from == 0;

				if (bulge)
				{
					if (bridge.type == parallel)
						bulge = ((jb1 - bridge.je < 6 and ib1 - bridge.ie < 3) or jb1 - bridge.je < 3) and
							residues[bridge.je].chain == residues[jb1].chain;
					else
						bulge = ((bridge.jb - je1 < 6 and ib1 - bridge.ie < 3) or bridge.jb - je1 < 3) and
							residues[je1].chain == residues[bridge.jb].chain;
				}
				
				if (bulge)
				{
					bridge.to = j;
					bridges[j].from = i;
				}
			}
		}
		
		for (uint32 i = 0; i < bridges.size(); ++i)
		{
			MBridge& bridge = bridges[i];
			if (bridge.from == 0)
			{
				for (uint32 j = i; j != 0; j = bridges[j].to)
					bridge.link.insert(j);
				for (uint32 j = bridge.to; j != 0; j = bridges[j].to)
					bridges[j].link = bridge.link;
			}
		}
		
		// Sheet
		
		set<uint16> sheetset, ladderset;
		for (uint32 i = 0; i < bridges.size(); ++i)
			ladderset.insert(i);
		
		char sheetname = 'A' - 1;
		char laddername = 'a';
		
		while (not ladderset.empty())
		{
			++sheetname;
			
			sheetset = bridges[*ladderset.begin()].link;
			bool done = false;
			while (not done)
			{
				done = true;
				for (uint32 l1 = 1; l1 < bridges.size(); ++l1)
				{
					foreach (uint32 l2, ladderset)
					{
						if ((bridges[l1].ie >= bridges[l2].ib and bridges[l1].ib <= bridges[l2].ie) or
							(bridges[l1].ie >= bridges[l2].jb and bridges[l1].ib <= bridges[l2].je) or
							(bridges[l1].je >= bridges[l2].ib and bridges[l1].jb <= bridges[l2].ie) or
							(bridges[l1].je >= bridges[l2].jb and bridges[l1].jb <= bridges[l2].je))
						{
							sheetset.insert(bridges[l2].link.begin(), bridges[l2].link.end());
							foreach (uint16 l, bridges[l2].link)
								ladderset.erase(l);
							done = false;
						}
					}
				}
			}
			
			for (uint32 i = 0; i < bridges.size(); ++i)
			{
				MBridge& bridge = bridges[i];
				if (sheetset.count(i) and bridge.from == 0)
				{
					if (++laddername > 'z')
						laddername = 'a';
					if (bridge.type == parallel)
						bridge.ladder = laddername;
					else
						bridge.ladder = toupper(laddername);
					bridge.sheet = sheetname;
					bridge.link = sheetset;
					for (uint32 j = bridge.to; j != 0; j = bridges[j].to)
					{
						bridges[j].ladder = laddername;
						bridges[j].sheet = sheetname;
						bridges[j].link = sheetset;
					}
				}
			}
		}
		
//		MResidue::Sheet(bridges);
//		MResidue::Markstrands(bridges);
	}
}

void MProtein::Center()
{
	vector<point> p;
	GetPoints(p);
	
	point t = center_points(p);
	
	Translate(point(-t.m_x, -t.m_y, -t.m_z));
}

MResidue& MProtein::GetResidue(char inChainID, uint16 inSeqNumber)
{
	MChain& chain = mChains[inChainID];
	if (chain.mResidues.empty())
		throw mas_exception(boost::format("Invalid chain id '%c'") % inChainID);
	return chain.GetResidueBySeqNumber(inSeqNumber);
}

void MProtein::GetCAlphaLocations(char inChain, vector<point>& outPoints) const
{
	if (inChain == 0)
		inChain = mChains.begin()->first;
	
	map<char,MChain>::const_iterator chain = mChains.find(inChain);
	assert(chain != mChains.end());
	
	foreach (const MResidue* r, chain->second.mResidues)
		outPoints.push_back(r->GetCAlpha());
}

point MProtein::GetCAlphaPosition(char inChain, int16 inPDBResSeq) const
{
	if (inChain == 0)
		inChain = mChains.begin()->first;
	
	map<char,MChain>::const_iterator chain = mChains.find(inChain);
	assert(chain != mChains.end());
	
	point result;
	
	foreach (const MResidue* r, chain->second.mResidues)
	{
		if (r->mSeqNumber != inPDBResSeq)
			continue;
		
		result = r->GetCAlpha();
	}
	
	return result;
}

void MProtein::CalculateSSBridges()
{
	
}

void MProtein::GetSequence(char inChain, entry& outEntry) const
{
	if (inChain == 0)
		inChain = mChains.begin()->first;
	
	map<char,MChain>::const_iterator chain = mChains.find(inChain);
	assert(chain != mChains.end());
	
	string seq;
	
	foreach (const MResidue* r, chain->second.mResidues)
	{
		foreach (const MAtom& a, r->mAtoms)
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
	
	foreach (const MResidue* r, chain->second.mResidues)
	{
		foreach (const MAtom& a, r->mAtoms)
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

void MProtein::WriteDSSP(ostream& os)
{
	boost::format kDSSPResidueLine(
	"%5.5d        !*             0   0    0      0, 0.0     0, 0.0     0, 0.0     0, 0.0   0.000 360.0 360.0 360.0 360.0    0.0    0.0    0.0");  

	for (map<char,MChain>::iterator c = mChains.begin(); c != mChains.end(); ++c)
	{
		c->second.WriteDSSP(os);
		if (next(c) != mChains.end())
			os << (kDSSPResidueLine % (c->second.mResidues.back()->mNumber + 1)) << endl;
	}
}

// --------------------------------------------------------------------

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
		throw logic_error("Protein A and B should have the same number of c-alpha atoms");
	
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

void create_entries(MProtein& a, MProtein& b, char chainA, char chainB,
	entry& ea, entry& eb)
{
	a.GetSequence(chainA, ea.m_seq);
	b.GetSequence(chainB, eb.m_seq);
	
	ea.m_positions = vector<int16>(ea.m_seq.length());
	eb.m_positions = vector<int16>(eb.m_seq.length());

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

	int32 ai = highA, bi = highB;
	int16 posNr = max(N, M) + 1;
	
	while (ai > 0 and bi > 0)
	{
		switch (traceback(ai, bi))
		{
			case 0:
				if (Distance(cAlphaA[ai], cAlphaB[bi]) < 3.5)
				{
					assert(posNr > 0);
					ea.m_positions[ai] = posNr;
					eb.m_positions[bi] = posNr;
					--posNr;
				}
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
}

void align_structures(const string& structureA, const string& structureB,
	uint32 iterations,
	substitution_matrix_family& mat, float gop, float gep, float magic)
{
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
	MProtein a(file_a, false);
	
	if (chainA == 0)
		chainA = a.mChains.begin()->first;

	stringstream file_b(pdb->GetDocument(pdbid_b));
	MProtein b(file_b, false);

	if (chainB == 0)
		chainB = b.mChains.begin()->first;

	align_proteins(a, chainA, b, chainB, iterations, mat, gop, gep, magic);
	
	if (VERBOSE)
		cerr << "RMSD: " << CalculateRMSD(a, b, chainA, chainB) << endl;

	entry ea(1, a.mID), eb(2, b.mID);
	create_entries(a, b, chainA, chainB, ea, eb);
	
	vector<entry*> aa, ab, ac;
	aa.push_back(&ea);
	ab.push_back(&eb);
	joined_node n(new leaf_node(ea), new leaf_node(eb), 0.1, 0.1);

	align(&n, aa, ab, ac, mat, gop, gep, magic, false);
	report(ac, cout, "clustalw");
	
	MProtein c;

	c.mChains['A'] = a.mChains[chainA];
	c.mChains['A'].SetChainID('A');

	c.mChains['B'] = b.mChains[chainB];
	c.mChains['B'].SetChainID('B');
	
	ofstream file_o(pdbid_a + chainA + '-' + pdbid_b + chainB + ".pdb");
	c.WritePDB(file_o);
}

// --------------------------------------------------------------------

void test_ss(const string& inID)
{
	CDatabankPtr pdb = LoadDatabank("pdb");
	
	stringstream file(pdb->GetDocument(inID));
	
	MProtein a(file, false);

	a.CalculateSecondaryStructure();
	
//	a.WritePDB(cout);
	a.WriteDSSP(cout);
}
