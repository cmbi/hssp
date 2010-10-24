// structure related stuff

#include "mas.h"

#include <set>
#include <numeric>
#include <functional>

#include <boost/format.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/foreach.hpp>
#define foreach BOOST_FOREACH
#include <boost/algorithm/string.hpp>

#include "utils.h"
#include "structure.h"

using namespace std;
namespace ba = boost::algorithm;
//namespace bm = boost::math;

// --------------------------------------------------------------------

const double
	kSSBridgeDistance = 3.0,
	kMinimalDistance = 0.5,
	kMinimalCADistance = 9.0,
	kMinHBondEnergy = -9.9,
	kMaxHBondEnergy = -0.5,
	kCouplingConstant = -332 * 0.42 * 0.2,
	kMaxPeptideBondLength = 2.5;

const double
	kRadiusN = 1.65,
	kRadiusCA = 1.87,
	kRadiusC = 1.76,
	kRadiusO = 1.4,
	kRadiusSideAtom = 1.8,
	kRadiusWater = 1.4;

// --------------------------------------------------------------------

MAtomType MapElement(string inElement)
{
	ba::trim(inElement);
	
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

MResidueType MapResidue(string inName)
{
	ba::trim(inName);

	MResidueType result = kUnknownResidue;
	
	for (uint32 i = 0; i < kResidueTypeCount; ++i)
	{
		if (inName == kResidueInfo[i].name)
		{
			result = kResidueInfo[i].type;
			break;
		}
	}
	
	return result;
}

// --------------------------------------------------------------------

inline double ParseFloat(const string& s)
{
	return boost::lexical_cast<double>(ba::trim_copy(s));
}

void MAtom::WritePDB(ostream& os) const
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
		   mLoc.mX % mLoc.mY % mLoc.mZ % mOccupancy % mTempFactor % mElement % charge) << endl;
}

const MResidueInfo kResidueInfo[] = {
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

struct MBridge
{
	MBridgeType		type;
	uint32			sheet, ladder;
	set<MBridge*>	link;
	deque<uint32>	i, j;
	char			chainI, chainJ;
	
	bool			operator<(const MBridge& b) const		{ return chainI < b.chainI or (chainI == b.chainI and i.front() < b.i.front()); }
};

// return true if any of the residues in bridge a is identical to any of the residues in bridge b
bool Linked(const MBridge& a, const MBridge& b)
{
	bool result = find_first_of(a.i.begin(), a.i.end(), b.i.begin(), b.i.end()) != a.i.end();
	if (result == false)
		result = find_first_of(a.i.begin(), a.i.end(), b.j.begin(), b.j.end()) != a.i.end();
	if (result == false)
		result = find_first_of(a.j.begin(), a.j.end(), b.i.begin(), b.i.end()) != a.j.end();
	if (result == false)
		result = find_first_of(a.j.begin(), a.j.end(), b.j.begin(), b.j.end()) != a.j.end();
	return result;
}

// --------------------------------------------------------------------

MResidue::MResidue(MChain& chain, uint32 inNumber,
		MResidue* inPrevious, const vector<MAtom>& inAtoms)
	: mChain(chain)
	, mPrev(inPrevious)
	, mNext(nil)
	, mSeqNumber(inAtoms.front().mResSeq)
	, mNumber(inNumber)
	, mType(MapResidue(inAtoms.front().mResName))
	, mSSBridgeNr(0)
	, mAccessibility(0)
	, mSecondaryStructure(loop)
	, mSheet(0)
{
	if (mPrev != nil)
		mPrev->mNext = this;
	
	fill(mHelixFlags, mHelixFlags + 3, helixNone);
	
	mBetaPartner[0].residue = mBetaPartner[1].residue = nil;
	
	mHBondDonor[0].energy = mHBondDonor[1].energy = mHBondAcceptor[0].energy = mHBondAcceptor[1].energy = 0;
	mHBondDonor[0].residue = mHBondDonor[1].residue = mHBondAcceptor[0].residue = mHBondAcceptor[1].residue = nil;

	static const MAtom kNullAtom = {};
	mN = mCA = mC = mO = kNullAtom;
	
	foreach (const MAtom& atom, inAtoms)
	{
		if (MapResidue(atom.mResName) != mType)
			throw mas_exception("inconsistent residue types in atom records");
		
		if (atom.mResSeq != mSeqNumber)
			throw mas_exception("inconsistent residue sequence numbers");
		
		if (atom.GetName() == " N  ")
			mN = atom;
		else if (atom.GetName() == " CA ")
			mCA = atom;
		else if (atom.GetName() == " C  ")
			mC = atom;
		else if (atom.GetName() == " O  ")
			mO = atom;
		else
			mSideChain.push_back(atom);
	}
	
	// assign the Hydrogen
	mH = GetN();
	
	if (mType != kProline and mPrev != nil)
	{
		const MAtom& pc = mPrev->GetC();
		const MAtom& po = mPrev->GetO();
		
		double CODistance = Distance(pc, po);
		
		mH.mLoc.mX += (pc.mLoc.mX - po.mLoc.mX) / CODistance; 
		mH.mLoc.mY += (pc.mLoc.mY - po.mLoc.mY) / CODistance; 
		mH.mLoc.mZ += (pc.mLoc.mZ - po.mLoc.mZ) / CODistance; 
	}
	
	// update the box containing all atoms
	mBox[0].mX = mBox[0].mY = mBox[0].mZ =  numeric_limits<double>::max();
	mBox[1].mX = mBox[1].mY = mBox[1].mZ = -numeric_limits<double>::max();
	
	ExtendBox(mN, kRadiusN + 2 * kRadiusWater);
	ExtendBox(mCA, kRadiusCA + 2 * kRadiusWater);
	ExtendBox(mC, kRadiusC + 2 * kRadiusWater);
	ExtendBox(mO, kRadiusO + 2 * kRadiusWater);
	foreach (const MAtom& atom, mSideChain)
		ExtendBox(atom, kRadiusSideAtom + 2 * kRadiusWater);
}

void MResidue::SetChainID(char inID)
{
	mC.SetChainID(inID);
	mCA.SetChainID(inID);
	mO.SetChainID(inID);
	mN.SetChainID(inID);
	mH.SetChainID(inID);
	for_each(mSideChain.begin(), mSideChain.end(), boost::bind(&MAtom::SetChainID, _1, inID));
}

char MResidue::GetChainID() const
{
	return mChain.GetChainID();
}

bool MResidue::ValidDistance(const MResidue& inNext) const
{
	return Distance(GetC(), inNext.GetN()) <= kMaxPeptideBondLength;
}

bool MResidue::TestBond(const MResidue* other) const
{
	return
		(mHBondAcceptor[0].residue == other and mHBondAcceptor[0].energy < kMaxHBondEnergy) or
		(mHBondAcceptor[1].residue == other and mHBondAcceptor[1].energy < kMaxHBondEnergy);
}

double MResidue::Phi() const
{
	double result = 360;
	if (mPrev != nil)
		result = DihedralAngle(mPrev->GetC(), GetN(), GetCAlpha(), GetC());
	return result;
}

double MResidue::Psi() const
{
	double result = 360;
	if (mNext != nil)
		result = DihedralAngle(GetN(), GetCAlpha(), GetC(), mNext->GetN());
	return result;
}

tr1::tuple<double,char> MResidue::Alpha() const
{
	double alhpa = 360;
	char chirality = ' ';
	
	const MResidue* nextNext = mNext ? mNext->Next() : nil;
	if (mPrev != nil and nextNext != nil)
	{
		alhpa = DihedralAngle(mPrev->GetCAlpha(), GetCAlpha(), mNext->GetCAlpha(), nextNext->GetCAlpha());
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
	const MResidue* prevPrev = mPrev ? mPrev->Prev() : nil;
	const MResidue* nextNext = mNext ? mNext->Next() : nil;
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
	if (mPrev != nil)
		result = CosinusAngle(GetC(), GetO(), mPrev->GetC(), mPrev->GetO());
	return result;
}

void MResidue::SetBetaPartner(uint32 n,
	MResidue* inResidue, uint32 inLadder, bool inParallel)
{
	assert(n == 0 or n == 1);
	
	mBetaPartner[n].residue = inResidue;
	mBetaPartner[n].ladder = inLadder;
	mBetaPartner[n].parallel = inParallel;
}

MBridgeParner MResidue::GetBetaPartner(uint32 n) const
{
	assert(n == 0 or n == 1);
	return mBetaPartner[n];
}

MHelixFlag MResidue::GetHelixFlag(uint32 inHelixStride) const
{
	assert(inHelixStride == 3 or inHelixStride == 4 or inHelixStride == 5);
	return mHelixFlags[inHelixStride - 3];
}

bool MResidue::IsHelixStart(uint32 inHelixStride) const
{
	assert(inHelixStride == 3 or inHelixStride == 4 or inHelixStride == 5);
	return mHelixFlags[inHelixStride - 3] == helixStart or mHelixFlags[inHelixStride - 3] == helixStartAndEnd;
}

void MResidue::SetHelixFlag(uint32 inHelixStride, MHelixFlag inHelixFlag)
{
	assert(inHelixStride == 3 or inHelixStride == 4 or inHelixStride == 5);
	mHelixFlags[inHelixStride - 3] = inHelixFlag;
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
{										// I.	a	d	II.	a	d		parallel    
	const MResidue* a = mPrev;			//		  \			  /
	const MResidue* b = this;			//		b	e		b	e
	const MResidue* c = mNext;			// 		  /			  \                      ..
	const MResidue* d = test->mPrev;	//		c	f		c	f
	const MResidue* e = test;			//
	const MResidue* f = test->mNext;	// III.	a <- f	IV. a	  f		antiparallel
										//		                                   
	MBridgeType result = btNoBridge;	//		b	 e      b <-> e                  
	if (a and c and d and f)			//                                          
	{									//		c -> d		c     d
		if ((TestBond(c, e) and TestBond(e, a)) or (TestBond(f, b) and TestBond(b, d)))
			result = btParallel;
		else if ((TestBond(c, d) and TestBond(f, a)) or (TestBond(e, b) and TestBond(b, e)))
			result = btAntiParallel;
	}
	
	return result;
}

void MResidue::ExtendBox(const MAtom& atom, double inRadius)
{
	if (mBox[0].mX > atom.mLoc.mX - inRadius)
		mBox[0].mX = atom.mLoc.mX - inRadius;
	if (mBox[0].mY > atom.mLoc.mY - inRadius)
		mBox[0].mY = atom.mLoc.mY - inRadius;
	if (mBox[0].mZ > atom.mLoc.mZ - inRadius)
		mBox[0].mZ = atom.mLoc.mZ - inRadius;
	if (mBox[1].mX < atom.mLoc.mX + inRadius)
		mBox[1].mX = atom.mLoc.mX + inRadius;
	if (mBox[1].mY < atom.mLoc.mY + inRadius)
		mBox[1].mY = atom.mLoc.mY + inRadius;
	if (mBox[1].mZ < atom.mLoc.mZ + inRadius)
		mBox[1].mZ = atom.mLoc.mZ + inRadius;
}

inline
bool MResidue::AtomIntersectsBox(const MAtom& atom, double inRadius) const
{
	return
		atom.mLoc.mX + inRadius >= mBox[0].mX and atom.mLoc.mX - inRadius <= mBox[1].mX and
		atom.mLoc.mY + inRadius >= mBox[0].mY and atom.mLoc.mY - inRadius <= mBox[1].mY and
		atom.mLoc.mZ + inRadius >= mBox[0].mZ and atom.mLoc.mZ - inRadius <= mBox[1].mZ;	
}

double MResidue::CalculateSurface(const vector<MPoint>& inPolyeder,
	const vector<double>& inWeights, const vector<MResidue*>& inResidues)
{
	double surface = CalculateSurface(mN, kRadiusN, inPolyeder, inWeights, inResidues) +
					 CalculateSurface(mCA, kRadiusCA, inPolyeder, inWeights, inResidues) +
					 CalculateSurface(mC, kRadiusC, inPolyeder, inWeights, inResidues) +
					 CalculateSurface(mO, kRadiusO, inPolyeder, inWeights, inResidues);
	
	foreach (const MAtom& atom, mSideChain)
		surface += CalculateSurface(atom, kRadiusSideAtom, inPolyeder, inWeights, inResidues);

	mAccessibility = floor(surface + 0.5);
	
	return surface;
}

class MAccumulator
{
  public:

	struct candidate
	{
		MPoint	location;
		double	radius;
		double	distance;
		
		bool operator<(const candidate& rhs) const
				{ return distance < rhs.distance; }
	};
	
	void operator()(const MPoint& a, const MPoint& b, double d, double r)
	{
		double distance = DistanceSquared(a, b);
		
		d += kRadiusWater;
		r += kRadiusWater;
		
		double test = d + r;
		test *= test;

		if (distance < test and distance > 0.0001)
		{
			candidate c = { b - a, r * r, distance };
			
			m_x.push_back(c);
			push_heap(m_x.begin(), m_x.end());
		}
	}
	
	void sort()
	{
		sort_heap(m_x.begin(), m_x.end());
	}

	vector<candidate>	m_x;
};

double MResidue::CalculateSurface(const MAtom& inAtom, double inRadius,
	const vector<MPoint>& inPolyeder, const vector<double>& inWeights, const vector<MResidue*>& inResidues)
{
	MAccumulator accumulate;
	
	foreach (MResidue* r, inResidues)
	{
		if (r->AtomIntersectsBox(inAtom, inRadius))
//		if (Distance(r->mCA, inAtom) < 10.0 + 2 * kRadiusWater)
		{
//cerr << "add residue " << r->mNumber << endl;
			accumulate(inAtom, r->mN, inRadius, kRadiusN);
			accumulate(inAtom, r->mCA, inRadius, kRadiusCA);
			accumulate(inAtom, r->mC, inRadius, kRadiusC);
			accumulate(inAtom, r->mO, inRadius, kRadiusO);
				
			foreach (const MAtom& atom, r->mSideChain)
				accumulate(inAtom, atom, inRadius, kRadiusSideAtom);
		}
	}

	accumulate.sort();

	double radius = inRadius + kRadiusWater;
	double surface = 0;
	
	for (uint32 i = 0; i < inPolyeder.size(); ++i)
	{
		MPoint xx = inPolyeder[i] * radius;
		
		bool free = true;
		for (uint32 k = 0; free and k < accumulate.m_x.size(); ++k)
			free = DistanceSquared(xx, accumulate.m_x[k].location) >= accumulate.m_x[k].radius;
		
		if (free)
			surface += inWeights[i];
	}
	
	return surface * radius * radius;
}

void MResidue::Translate(const MPoint& inTranslation)
{
	mN.Translate(inTranslation);
	mCA.Translate(inTranslation);
	mC.Translate(inTranslation);
	mO.Translate(inTranslation);
	mH.Translate(inTranslation);
	for_each(mSideChain.begin(), mSideChain.end(), boost::bind(&MAtom::Translate, _1, inTranslation));
}

void MResidue::Rotate(const MQuaternion& inRotation)
{
	mN.Rotate(inRotation);
	mCA.Rotate(inRotation);
	mC.Rotate(inRotation);
	mO.Rotate(inRotation);
	mH.Rotate(inRotation);
	for_each(mSideChain.begin(), mSideChain.end(), boost::bind(&MAtom::Rotate, _1, inRotation));
}

void MResidue::GetPoints(vector<MPoint>& outPoints) const
{
	outPoints.push_back(mN);
	outPoints.push_back(mCA);
	outPoints.push_back(mC);
	outPoints.push_back(mO);
	foreach (const MAtom& a, mSideChain)
		outPoints.push_back(a);
}

void MResidue::WritePDB(std::ostream& os)
{
	mN.WritePDB(os);
	mCA.WritePDB(os);
	mC.WritePDB(os);
	mO.WritePDB(os);
	
	for_each(mSideChain.begin(), mSideChain.end(), boost::bind(&MAtom::WritePDB, _1, ref(os)));
}

// --------------------------------------------------------------------

void MChain::SetChainID(char inID)
{
	mChainID = inID;
	for_each(mResidues.begin(), mResidues.end(), boost::bind(&MResidue::SetChainID, _1, inID));
}

void MChain::Translate(const MPoint& inTranslation)
{
	for_each(mResidues.begin(), mResidues.end(), boost::bind(&MResidue::Translate, _1, inTranslation));
}

void MChain::Rotate(const MQuaternion& inRotation)
{
	for_each(mResidues.begin(), mResidues.end(), boost::bind(&MResidue::Rotate, _1, inRotation));
}

void MChain::WritePDB(std::ostream& os)
{
	for_each(mResidues.begin(), mResidues.end(), boost::bind(&MResidue::WritePDB, _1, ref(os)));
	
	boost::format ter("TER    %4.4d      %3.3s %c%4.4d%c");
	
	MResidue* last = mResidues.back();
	
	os << (ter % (last->GetCAlpha().mSerial + 1) % kResidueInfo[last->GetType()].name % mChainID % last->GetNumber() % ' ') << endl;
}

MResidue& MChain::GetResidueBySeqNumber(uint16 inSeqNumber)
{
	vector<MResidue*>::iterator r = find_if(mResidues.begin(), mResidues.end(),
		boost::bind(&MResidue::GetSeqNumber, _1) == inSeqNumber);
	if (r == mResidues.end())
		throw mas_exception(boost::format("Residue %d not found") % inSeqNumber);
	return **r;
}

// --------------------------------------------------------------------

MProtein::MProtein(istream& is, bool cAlphaOnly)
	: mResidueCount(0)
{
	bool model = false;
	vector<MAtom> atoms;
	
	while (not is.eof())
	{
		string line;
		getline(is, line);
		
		if (ba::starts_with(line, "HEADER"))
		{
			mHeader = line.substr(0, 62);
			ba::trim(mHeader);
			mID = line.substr(62, 4);
			continue;
		}

		if (ba::starts_with(line, "COMPND"))
		{
			mCompound.push_back(ba::trim_copy(line));
			continue;
		}

		if (ba::starts_with(line, "SOURCE"))
		{
			mSource.push_back(ba::trim_copy(line));
			continue;
		}
		
		if (ba::starts_with(line, "AUTHOR"))
		{
			mAuthor.push_back(ba::trim_copy(line));
			continue;
		}
		
		// brain dead support for only the first model in the file (NMR)
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
		
		if (ba::starts_with(line, "TER   "))
		{
			if (atoms.empty())
				throw mas_exception("no atoms read before TER record");
			
			AddResidue(atoms);
			atoms.clear();
			continue;
		}
		
		if (ba::starts_with(line, "ATOM  ") or ba::starts_with(line, "HETATM"))
			//	1 - 6	Record name "ATOM "
		{
			if (cAlphaOnly and line.substr(12, 4) != " CA ")
				continue;

			MAtom atom = {};

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
			atom.mLoc.mX = ParseFloat(line.substr(30, 8));
			//	39 - 46	Real(8.3) y Orthogonal coordinates for Y in Angstroms.
			atom.mLoc.mY = ParseFloat(line.substr(38, 8));
			//	47 - 54	Real(8.3) z Orthogonal coordinates for Z in Angstroms.
			atom.mLoc.mZ = ParseFloat(line.substr(46, 8));
			//	55 - 60	Real(6.2) occupancy Occupancy.
			atom.mOccupancy = ParseFloat(line.substr(54, 6));
			//	61 - 66	Real(6.2) tempFactor Temperature factor.
			atom.mTempFactor = ParseFloat(line.substr(60, 6));
			//	77 - 78	LString(2) element Element symbol, right-justified.
			if (line.length() > 76)
				line.copy(atom.mElement, 2, 76);
			//	79 - 80	LString(2) charge Charge on the atom.
			atom.mCharge = 0;
			
			try
			{
				atom.mType = MapElement(line.substr(77, 2));
			}
			catch (exception& e)
			{
				if (VERBOSE)
					cerr << e.what() << endl;
				atom.mType = kUnknownAtom;
			}
			
			if (atom.mAltLoc != ' ' and atom.mAltLoc != 'A')
			{
				if (VERBOSE)
					cerr << "skipping alternate atom record " << atom.mResName << endl;
				continue;
			}
			
			if (not atoms.empty() and
				(atom.mResSeq != atoms.back().mResSeq or (atom.mResSeq == atoms.back().mResSeq and atom.mICode != atoms.back().mICode)))
			{
				AddResidue(atoms);
				atoms.clear();
			}

			if (atom.mType != kHydrogen)
				atoms.push_back(atom);
		}
	}
	
	mChains.erase(
		remove_if(mChains.begin(), mChains.end(), boost::bind(&MChain::Empty, _1)),
		mChains.end());
}

string MProtein::GetCompound() const
{
	string result = "(missing)";
	if (not mCompound.empty())
	{
		result = mCompound.front();
		if (ba::starts_with(result.substr(10), "MOL_ID: ") and mCompound.size() > 1)
			result = mCompound[1];
	}
	return result;
}

string MProtein::GetSource() const
{
	string result = "(missing)";
	if (not mSource.empty())
	{
		result = mSource.front();
		if (ba::starts_with(result.substr(10), "MOL_ID: ") and mSource.size() > 1)
			result = mSource[1];
	}
	return result;
}

string MProtein::GetAuthor() const
{
	string result = "(missing)";
	if (not mAuthor.empty())
	{
		result = mAuthor.front();
		if (ba::starts_with(result.substr(10), "MOL_ID: ") and mAuthor.size() > 1)
			result = mAuthor[1];
	}
	return result;
}

void MProtein::GetStatistics(uint32& outNrOfResidues, uint32& outNrOfChains,
	uint32& outNrOfSSBridges, uint32& outNrOfIntraChainSSBridges) const
{
	outNrOfResidues = mResidueCount;
	outNrOfChains = mChains.size();
	outNrOfSSBridges = mSSBonds.size();
	
	outNrOfIntraChainSSBridges = 0;
	for (vector<pair<MResidueID,MResidueID>>::const_iterator ri = mSSBonds.begin(); ri != mSSBonds.end(); ++ri)
	{
		if (ri->first.chain == ri->second.chain)
			++outNrOfIntraChainSSBridges;
	}
}

void MProtein::AddResidue(const vector<MAtom>& inAtoms)
{
	MChain& chain = GetChain(inAtoms.front().mChainID);
	vector<MResidue*>& residues(chain.GetResidues());

	MResidue* prev = nil;
	if (not residues.empty())
		prev = residues.back();

	bool hasN = false, hasCA = false, hasC = false, hasO = false;
	foreach (const MAtom& atom, inAtoms)
	{
		if (atom.GetName() == " N  ")
			hasN = true;
		if (atom.GetName() == " CA ")
			hasCA = true;
		if (atom.GetName() == " C  ")
			hasC = true;
		if (atom.GetName() == " O  ")
			hasO = true;
	}
	
	if (hasN and hasCA and hasC and hasO)
	{
		uint32 resNumber = mResidueCount + mChains.size();
		residues.push_back(new MResidue(chain, resNumber, prev, inAtoms));
		++mResidueCount;
	}
	else if (VERBOSE)
		cerr << "ignoring incomplete residue " << inAtoms.front().mResName << " (" << inAtoms.front().mResSeq << ')' << endl;
}

const MChain& MProtein::GetChain(char inChainID) const
{
	for (uint32 i = 0; i < mChains.size(); ++i)
		if (mChains[i]->GetChainID() == inChainID)
			return *mChains[i];
	
	throw mas_exception("Chain not found");
	return *mChains.front();
}

MChain& MProtein::GetChain(char inChainID)
{
	for (uint32 i = 0; i < mChains.size(); ++i)
		if (mChains[i]->GetChainID() == inChainID)
			return *mChains[i];
	
	mChains.push_back(new MChain(inChainID));
	return *mChains.back();
}

void MProtein::GetPoints(std::vector<MPoint>& outPoints) const
{
	foreach (const MChain* chain, mChains)
	{
		foreach (const MResidue* r, chain->GetResidues())
			r->GetPoints(outPoints);
	}
}

void MProtein::Translate(const MPoint& inTranslation)
{
	foreach (MChain* chain, mChains)
		chain->Translate(inTranslation);
}

void MProtein::Rotate(const MQuaternion& inRotation)
{
	foreach (MChain* chain, mChains)
		chain->Rotate(inRotation);
}

void MProtein::CalculateSecondaryStructure()
{
	if (VERBOSE)
		cerr << "Assigning sulfur bridges" << endl;
	
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

	vector<MResidue*> residues;
	residues.reserve(mResidueCount);
	foreach (const MChain* chain, mChains)
		residues.insert(residues.end(), chain->GetResidues().begin(), chain->GetResidues().end());
	
	if (VERBOSE)
		cerr << "using " << residues.size() << " residues" << endl;

	CalculateHBondEnergies(residues);
	CalculateBetaSheets(residues);
	CalculateAlphaHelices(residues);
	CalculateAccessibilities(residues);
}

void MProtein::CalculateHBondEnergies(std::vector<MResidue*> inResidues)
{
	if (VERBOSE)
		cerr << "Calculate H-bond energies" << endl;
	
	// Calculate the HBond energies
	for (uint32 i = 0; i < inResidues.size() - 1; ++i)
	{
		MResidue* ri = inResidues[i];
		
		for (uint32 j = i + 1; j < inResidues.size(); ++j)
		{
			MResidue* rj = inResidues[j];
			
			if (Distance(ri->GetCAlpha(), rj->GetCAlpha()) < kMinimalCADistance)
			{
				MResidue::CalculateHBondEnergy(*ri, *rj);
				if (j != i + 1)
					MResidue::CalculateHBondEnergy(*rj, *ri);
			}
		}
	}
}

void MProtein::CalculateAlphaHelices(std::vector<MResidue*> inResidues)
{
	if (VERBOSE)
		cerr << "Calculate alhpa helices" << endl;
	
	// Helix and Turn
	foreach (const MChain* chain, mChains)
	{
		for (uint32 stride = 3; stride <= 5; ++stride)
		{
			vector<MResidue*> res(chain->GetResidues());
			if (res.size() < stride)
				continue;
			
			for (uint32 i = 0; i < res.size() - stride; ++i)
			{
				if (MResidue::TestBond(res[i + stride], res[i]))
				{
					res[i + stride]->SetHelixFlag(stride, helixEnd);
					for (uint32 j = i + 1; j < i + stride; ++j)
					{
						if (res[j]->GetHelixFlag(stride) == helixNone)
							res[j]->SetHelixFlag(stride, helixMiddle);
					}
					
					if (res[i]->GetHelixFlag(stride) == helixEnd)
						res[i]->SetHelixFlag(stride, helixStartAndEnd);
					else
						res[i]->SetHelixFlag(stride, helixStart);
				}
			}
		}
	}
	
	foreach (MResidue* r, inResidues)
	{
		double kappa = r->Kappa();
		r->SetBend(kappa != 360 and kappa > 70);
	}

	for (uint32 i = 1; i < inResidues.size() - 4; ++i)
	{
		if (inResidues[i]->IsHelixStart(4) and inResidues[i - 1]->IsHelixStart(4))
		{
			for (uint32 j = i; j <= i + 3; ++j)
				inResidues[j]->SetSecondaryStructure(alphahelix);
		}
	}

	for (uint32 i = 1; i < inResidues.size() - 3; ++i)
	{
		if (inResidues[i]->IsHelixStart(3) and inResidues[i - 1]->IsHelixStart(3))
		{
			bool empty = true;
			for (uint32 j = i; empty and j <= i + 2; ++j)
				empty = inResidues[j]->GetSecondaryStructure() == loop or inResidues[j]->GetSecondaryStructure() == helix_3;
			if (empty)
			{
				for (uint32 j = i; j <= i + 2; ++j)
					inResidues[j]->SetSecondaryStructure(helix_3);
			}
		}
	}

	for (uint32 i = 1; i < inResidues.size() - 5; ++i)
	{
		if (inResidues[i]->IsHelixStart(5) and inResidues[i - 1]->IsHelixStart(5))
		{
			bool empty = true;
			for (uint32 j = i; empty and j <= i + 4; ++j)
				empty = inResidues[j]->GetSecondaryStructure() == loop or inResidues[j]->GetSecondaryStructure() == helix_5;
			if (empty)
			{
				for (uint32 j = i; j <= i + 4; ++j)
					inResidues[j]->SetSecondaryStructure(helix_5);
			}
		}
	}
			
	for (uint32 i = 1; i < inResidues.size() - 1; ++i)
	{
		if (inResidues[i]->GetSecondaryStructure() == loop)
		{
			bool isTurn = false;
			for (uint32 stride = 3; stride <= 5 and not isTurn; ++stride)
			{
				for (uint32 k = 1; k < stride and not isTurn; ++k)
					isTurn = (i >= k) and inResidues[i - k]->IsHelixStart(stride);
			}
			
			if (isTurn)
				inResidues[i]->SetSecondaryStructure(turn);
			else if (inResidues[i]->IsBend())
				inResidues[i]->SetSecondaryStructure(bend);
		}
	}
}

void MProtein::CalculateBetaSheets(std::vector<MResidue*> inResidues)
{
	if (VERBOSE)
		cerr << "Calculate beta sheets" << endl;
	
	// Calculate Bridges
	vector<MBridge> bridges;
	if (inResidues.size() > 4)
	{
		for (uint32 i = 1; i < inResidues.size() - 4; ++i)
		{
			MResidue* ri = inResidues[i];
			
			for (uint32 j = i + 3; j < inResidues.size() - 1; ++j)
			{
				MResidue* rj = inResidues[j];
				
				MBridgeType type = ri->TestBridge(rj);
				if (type == btNoBridge)
					continue;
				
				bool found = false;
				foreach (MBridge& bridge, bridges)
				{
					if (type != bridge.type or i != bridge.i.back() + 1)
						continue;
					
					if (type == btParallel and bridge.j.back() + 1 == j)
					{
						bridge.i.push_back(i);
						bridge.j.push_back(j);
						found = true;
						break;
					}
	
					if (type == btAntiParallel and bridge.j.front() - 1 == j)
					{
						bridge.i.push_back(i);
						bridge.j.push_front(j);
						found = true;
						break;
					}
				}
				
				if (not found)
				{
					MBridge bridge = {};
					
					bridge.type = type;
					bridge.i.push_back(i);
					bridge.chainI = ri->GetChainID();
					bridge.j.push_back(j);
					bridge.chainJ = rj->GetChainID();
					
					bridges.push_back(bridge);
				}
			}
		}
	}

	// extend ladders
	sort(bridges.begin(), bridges.end());
	
	for (uint32 i = 0; i < bridges.size(); ++i)
	{
		for (uint32 j = i + 1; j < bridges.size(); ++j)
		{
			if (bridges[i].type != bridges[j].type or
				bridges[i].chainI != bridges[j].chainI or
				bridges[i].chainJ != bridges[j].chainJ or
				bridges[j].i.front() - bridges[i].i.back() >= 6)
			{
				continue;
			}
			
//			uint32 ibi = bridges[i].i.front();
			uint32 iei = bridges[i].i.back();
			uint32 jbi = bridges[i].j.front();
			uint32 jei = bridges[i].j.back();
			uint32 ibj = bridges[j].i.front();
//			uint32 iej = bridges[j].i.back();
			uint32 jbj = bridges[j].j.front();
			uint32 jej = bridges[j].j.back();
			
			bool bulge;
			if (bridges[i].type == btParallel)
				bulge = (jbj - jei < 6 and ibj - iei < 3) or (jbj - jei < 3);
			else
				bulge = (jbi - jej < 6 and ibj - iei < 3) or (jbi - jej < 3);

			if (bulge)
			{
				bridges[i].i.insert(bridges[i].i.end(), bridges[j].i.begin(), bridges[j].i.end());
				if (bridges[i].type == btParallel)
					bridges[i].j.insert(bridges[i].j.end(), bridges[j].j.begin(), bridges[j].j.end());
				else
					bridges[i].j.insert(bridges[i].j.begin(), bridges[j].j.begin(), bridges[j].j.end());
				bridges.erase(bridges.begin() + j);
				--j;
			}
		}
	}

	// Sheet
	set<MBridge*> ladderset;
	foreach (MBridge& bridge, bridges)
		ladderset.insert(&bridge);
	
	uint32 sheet = 1, ladder = 0;
	while (not ladderset.empty())
	{
		set<MBridge*> sheetset;
		sheetset.insert(*ladderset.begin());
		ladderset.erase(ladderset.begin());

		bool done = false;
		while (not done)
		{
			done = true;
			foreach (MBridge* a, sheetset)
			{
				foreach (MBridge* b, ladderset)
				{
					if (Linked(*a, *b))
					{
						sheetset.insert(b);
						ladderset.erase(b);
						done = false;
						break;
					}
				}
				if (not done)
					break;
			}
		}

		foreach (MBridge* bridge, sheetset)
		{
			bridge->ladder = ladder;
			bridge->sheet = sheet;
			bridge->link = sheetset;
			
			++ladder;
		}
		
		++sheet;
	}

	foreach (MBridge& bridge, bridges)
	{
		// find out if any of the i and j set members already have
		// a bridge assigned, if so, we're assigning bridge 2
		
		uint32 betai = 0, betaj = 0;
		
		foreach (uint32 l, bridge.i)
		{
			if (inResidues[l]->GetBetaPartner(0).residue != nil)
			{
				betai = 1;
				break;
			}
		}

		foreach (uint32 l, bridge.j)
		{
			if (inResidues[l]->GetBetaPartner(0).residue != nil)
			{
				betaj = 1;
				break;
			}
		}
		
		MSecondaryStructure ss = betabridge;
		if (bridge.i.size() > 1)
			ss = strand;
		
		if (bridge.type == btParallel)
		{
			deque<uint32>::iterator j = bridge.j.begin();
			foreach (uint32 i, bridge.i)
				inResidues[i]->SetBetaPartner(betai, inResidues[*j++], bridge.ladder, true);

			j = bridge.i.begin();
			foreach (uint32 i, bridge.j)
				inResidues[i]->SetBetaPartner(betaj, inResidues[*j++], bridge.ladder, true);
		}
		else
		{
			deque<uint32>::reverse_iterator j = bridge.j.rbegin();
			foreach (uint32 i, bridge.i)
				inResidues[i]->SetBetaPartner(betai, inResidues[*j++], bridge.ladder, false);

			j = bridge.i.rbegin();
			foreach (uint32 i, bridge.j)
				inResidues[i]->SetBetaPartner(betaj, inResidues[*j++], bridge.ladder, false);
		}

		for (uint32 i = bridge.i.front(); i <= bridge.i.back(); ++i)
		{
			if (inResidues[i]->GetSecondaryStructure() != strand)
				inResidues[i]->SetSecondaryStructure(ss);
			inResidues[i]->SetSheet(bridge.sheet);
		}

		for (uint32 i = bridge.j.front(); i <= bridge.j.back(); ++i)
		{
			if (inResidues[i]->GetSecondaryStructure() != strand)
				inResidues[i]->SetSecondaryStructure(ss);
			inResidues[i]->SetSheet(bridge.sheet);
		}
	}
}

void CreateTriangle(const MPoint& p1, const MPoint& p2, const MPoint& p3, int level,
	vector<MPoint>& outPolyeders, vector<double>& outWeights)
{
	if (level > 0)
	{
		--level;
		MPoint p4 = p1 + p2;	p4.Normalize();
		MPoint p5 = p2 + p3;	p5.Normalize();
		MPoint p6 = p3 + p1;	p6.Normalize();
		
		CreateTriangle(p1, p4, p6, level, outPolyeders, outWeights);
		CreateTriangle(p4, p2, p5, level, outPolyeders, outWeights);
		CreateTriangle(p4, p5, p6, level, outPolyeders, outWeights);
		CreateTriangle(p5, p3, p6, level, outPolyeders, outWeights);
	}
	else
	{
		MPoint p = p1 + p2 + p3;
		p.Normalize();
		outPolyeders.push_back(p);
		
		p = CrossProduct(p3 - p1, p2 - p1);

		double l = p.Normalize();
		outWeights.push_back(l / 2);
	}
}

void MProtein::CalculateAccessibilities(std::vector<MResidue*> inResidues)
{
	if (VERBOSE)
		cerr << "Calculate accessibilities" << endl;

	// start by creating a icosahedron. Since this is constant, we can 
	// one day insert the constant data here.
	const double kXVertex = 0, kYVertex = 0.8506508, kZVertex = 0.5257311;
	
	const MPoint v[12] = {
		{ kXVertex, -kYVertex, -kZVertex }, { -kZVertex, kXVertex, -kYVertex }, { -kYVertex, -kZVertex, kXVertex },
		{ kXVertex, -kYVertex,  kZVertex }, {  kZVertex, kXVertex, -kYVertex }, { -kYVertex,  kZVertex, kXVertex },
		{ kXVertex,  kYVertex, -kZVertex }, { -kZVertex, kXVertex,  kYVertex }, {  kYVertex, -kZVertex, kXVertex },
		{ kXVertex,  kYVertex,  kZVertex }, {  kZVertex, kXVertex,  kYVertex }, {  kYVertex,  kZVertex, kXVertex },
	};
	
	const uint32 kOrder = 2;
	
	vector<MPoint> polyeder;
	vector<double> weights;
	
	for (uint32 i = 0; i < 10; ++i)
	{
		for (uint32 j = i + 1; j < 11; ++j)
		{
			if (Distance(v[i], v[j]) < 1.1)
			{
				for (uint32 k = j + 1; k < 12; ++k)
				{
					if (Distance(v[i], v[k]) < 1.1 and Distance(v[j], v[k]) < 1.1)
						CreateTriangle(v[i], v[j], v[k], kOrder, polyeder, weights);
				}
			}
		}
	}
	
	double a = std::accumulate(weights.begin(), weights.end(), 0.0);
	a = 4 * kPI / a;
	transform(weights.begin(), weights.end(), weights.begin(), boost::bind(multiplies<double>(), _1, a));

	foreach (MResidue* residue, inResidues)
		residue->CalculateSurface(polyeder, weights, inResidues);
}

void MProtein::Center()
{
	vector<MPoint> p;
	GetPoints(p);
	
	MPoint t = CenterPoints(p);
	
	Translate(MPoint(-t.mX, -t.mY, -t.mZ));
}

void MProtein::SetChain(char inChainID, const MChain& inChain)
{
	MChain& chain(GetChain(inChainID));
	chain = inChain;
	chain.SetChainID(inChainID);
}

MResidue& MProtein::GetResidue(char inChainID, uint16 inSeqNumber)
{
	MChain& chain = GetChain(inChainID);
	if (chain.GetResidues().empty())
		throw mas_exception(boost::format("Invalid chain id '%c'") % inChainID);
	return chain.GetResidueBySeqNumber(inSeqNumber);
}

void MProtein::GetCAlphaLocations(char inChain, vector<MPoint>& outPoints) const
{
	if (inChain == 0)
		inChain = mChains.front()->GetChainID();
	
	foreach (const MResidue* r, GetChain(inChain).GetResidues())
		outPoints.push_back(r->GetCAlpha());
}

MPoint MProtein::GetCAlphaPosition(char inChain, int16 inPDBResSeq) const
{
	if (inChain == 0)
		inChain = mChains.front()->GetChainID();
	
	MPoint result;
	foreach (const MResidue* r, GetChain(inChain).GetResidues())
	{
		if (r->GetSeqNumber() != inPDBResSeq)
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
		inChain = mChains.front()->GetChainID();
	
	string seq;
	foreach (const MResidue* r, GetChain(inChain).GetResidues())
	{
		seq += kResidueInfo[r->GetType()].code;
		outEntry.m_positions.push_back(r->GetSeqNumber());
	}
	
	outEntry.m_seq = encode(seq);
}

void MProtein::GetSequence(char inChain, sequence& outSequence) const
{
	if (inChain == 0)
		inChain = mChains.front()->GetChainID();
	
	string seq;
	foreach (const MResidue* r, GetChain(inChain).GetResidues())
		seq += kResidueInfo[r->GetType()].code;
	
	outSequence = encode(seq);
}

void MProtein::WritePDB(ostream& os)
{
	foreach (MChain* chain, mChains)
		chain->WritePDB(os);
}
