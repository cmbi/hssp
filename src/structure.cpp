// Copyright Maarten L. Hekkelman, Radboud University 2008-2011.
//   Distributed under the Boost Software License, Version 1.0.
//       (See accompanying file LICENSE_1_0.txt or copy at    
//             http://www.boost.org/LICENSE_1_0.txt)      
// 
// structure related stuff

#include "mas.h"

#include <set>
#include <numeric>
#include <functional>

#include <boost/bind.hpp>
#include <boost/format.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/foreach.hpp>
#define foreach BOOST_FOREACH
#include <boost/algorithm/string.hpp>
#include <boost/math/special_functions/round.hpp>

#include "align-2d.h"
#include "utils.h"
#include "buffer.h"
#include "structure.h"

using namespace std;
namespace ba = boost::algorithm;
namespace bm = boost::math;

// --------------------------------------------------------------------

const double
	kSSBridgeDistance = 3.0,
	kMinimalDistance = 0.5,
	kMinimalCADistance = 9.0,
	kMinHBondEnergy = -9.9,
	kMaxHBondEnergy = -0.5,
	kCouplingConstant = -27.888,	//	= -332 * 0.42 * 0.2
	kMaxPeptideBondLength = 2.5;

const double
	kRadiusN = 1.65,
	kRadiusCA = 1.87,
	kRadiusC = 1.76,
	kRadiusO = 1.4,
	kRadiusSideAtom = 1.8,
	kRadiusWater = 1.4;

// --------------------------------------------------------------------

namespace
{

// we use a fibonacci spheres to calculate the even distribution of the dots
class MSurfaceDots
{
  public:
	static MSurfaceDots&	Instance();
	
	uint32					size() const					{ return mPoints.size(); }
	const MPoint&			operator[](uint32 inIx) const	{ return mPoints[inIx]; }
	double					weight() const					{ return mWeight; }

  private:
							MSurfaceDots(int32 inN);

	vector<MPoint>			mPoints;
	double					mWeight;
};

MSurfaceDots& MSurfaceDots::Instance()
{
	const uint32 kN = 200;
	
	static MSurfaceDots sInstance(kN);
	return sInstance;
}

MSurfaceDots::MSurfaceDots(int32 N)
{
	int32 P = 2 * N + 1;
	
	const double
		kGoldenRatio = (1 + sqrt(5.0)) / 2;
	
	mWeight = (4 * kPI) / P;
	
	for (int32 i = -N; i <= N; ++i)
	{
		double lat = asin((2.0 * i) / P);
		double lon = fmod(i, kGoldenRatio) * 2 * kPI / kGoldenRatio;
		
		MPoint p;
		p.mX = sin(lon) * cos(lat);
		p.mY = cos(lon) * cos(lat);
		p.mZ = 			  sin(lat);

		mPoints.push_back(p);
	}
}

}

// --------------------------------------------------------------------

MAtomType MapElement(string inElement)
{
	ba::trim(inElement);
	ba::to_upper(inElement);
	
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
	else if (inElement == "CL")
		result = kChlorine;
	else if (inElement == "K")
		result = kPotassium;
	else if (inElement == "MG")
		result = kMagnesium;
	else if (inElement == "CA")
		result = kCalcium;
	else if (inElement == "ZN")
		result = kZinc;
	else if (inElement == "SE")
		result = kSelenium;
	else
		throw mas_exception(boost::format("Unsupported element '%1%'") % inElement);
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

MResidueType MapResidue(char inCode)
{
	MResidueType result = kUnknownResidue;
	
	for (uint32 i = 0; i < kResidueTypeCount; ++i)
	{
		if (inCode == kResidueInfo[i].code)
		{
			result = kResidueInfo[i].type;
			break;
		}
	}
	
	return result;
}

// --------------------------------------------------------------------
// a custom float parser, optimised for speed (and the way floats are represented in a PDB file)

double ParseFloat(const string& s)
{
	double result = 0;
	bool negate = false;
	double div = 10;
	
	enum State {
		pStart, pSign, pFirst, pSecond
	} state = pStart;
	
	for (string::const_iterator ch = s.begin(); ch != s.end(); ++ch)
	{
		switch (state)
		{
			case pStart:
				if (isspace(*ch))
					continue;
				if (*ch == '-')
				{
					negate = true;
					state = pSign;
				}
				else if (*ch == '+')
					state = pSign;
				else if (*ch == '.')
					state = pSecond;
				else if (isdigit(*ch))
				{
					result = *ch - '0';
					state = pFirst;
				}
				else
					throw mas_exception(boost::format("invalid formatted floating point number '%1%'") % s);
				break;
			
			case pSign:
				if (*ch == '.')
					state = pSecond;
				else if (isdigit(*ch))
				{
					state = pFirst;
					result = *ch - '0';
				}
				else
					throw mas_exception(boost::format("invalid formatted floating point number '%1%'") % s);
				break;
			
			case pFirst:
				if (*ch == '.')
					state = pSecond;
				else if (isdigit(*ch))
					result = 10 * result + (*ch - '0');
				else
					throw mas_exception(boost::format("invalid formatted floating point number '%1%'") % s);
				break;
			
			case pSecond:
				if (isdigit(*ch))
				{
					result += (*ch - '0') / div;
					div *= 10;
				}
				else
					throw mas_exception(boost::format("invalid formatted floating point number '%1%'") % s);
				break;
		}
	}
	
	if (negate)
		result = -result;
	
	return result;
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
	{ kUnknownResidue,	'X', "UNK" },
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

ostream& operator<<(ostream& os, const MBridge& b)
{
	os << '[' << (b.type == btParallel ? "p" : "a") << ':' << b.i.front() << '-' << b.i.back() << '/' << b.j.front() << '-' << b.j.back() << ']';
	return os;
}

// return true if any of the residues in bridge a is identical to any of the residues in bridge b
bool Linked(const MBridge& a, const MBridge& b)
{
	return
		find_first_of(a.i.begin(), a.i.end(), b.i.begin(), b.i.end()) != a.i.end() or
		find_first_of(a.i.begin(), a.i.end(), b.j.begin(), b.j.end()) != a.i.end() or
		find_first_of(a.j.begin(), a.j.end(), b.i.begin(), b.i.end()) != a.j.end() or
		find_first_of(a.j.begin(), a.j.end(), b.j.begin(), b.j.end()) != a.j.end();
}

// --------------------------------------------------------------------

MResidue::MResidue(uint32 inNumber,
		MResidue* inPrevious, const vector<MAtom>& inAtoms)
	: mChainID(0)
	, mPrev(inPrevious)
	, mNext(nullptr)
	, mSeqNumber(inAtoms.front().mResSeq)
	, mNumber(inNumber)
	, mInsertionCode(inAtoms.front().mICode)
	, mType(MapResidue(inAtoms.front().mResName))
	, mSSBridgeNr(0)
	, mAccessibility(0)
	, mSecondaryStructure(loop)
	, mSheet(0)
{
	if (mPrev != nullptr)
		mPrev->mNext = this;
	
	fill(mHelixFlags, mHelixFlags + 3, helixNone);
	
	mBetaPartner[0].residue = mBetaPartner[1].residue = nullptr;
	
	mHBondDonor[0].energy = mHBondDonor[1].energy = mHBondAcceptor[0].energy = mHBondAcceptor[1].energy = 0;
	mHBondDonor[0].residue = mHBondDonor[1].residue = mHBondAcceptor[0].residue = mHBondAcceptor[1].residue = nullptr;

	static const MAtom kNullAtom = {};
	mN = mCA = mC = mO = kNullAtom;
	
	foreach (const MAtom& atom, inAtoms)
	{
		if (mChainID == 0)
			mChainID = atom.mChainID;
		
		if (MapResidue(atom.mResName) != mType)
			throw mas_exception(
				boost::format("inconsistent residue types in atom records for residue %1% (%2% != %3%)")
					% inNumber % atom.mResName % inAtoms.front().mResName);
		
		if (atom.mResSeq != mSeqNumber)
			throw mas_exception(boost::format("inconsistent residue sequence numbers (%1% != %2%)") % atom.mResSeq % mSeqNumber);
		
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
	
	if (mType != kProline and mPrev != nullptr)
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
	
	mRadius = mBox[1].mX - mBox[0].mX;
	if (mRadius < mBox[1].mY - mBox[0].mY)
		mRadius = mBox[1].mY - mBox[0].mY;
	if (mRadius < mBox[1].mZ - mBox[0].mZ)
		mRadius = mBox[1].mZ - mBox[0].mZ;
	
	mCenter.mX = (mBox[0].mX + mBox[1].mX) / 2;
	mCenter.mY = (mBox[0].mY + mBox[1].mY) / 2;
	mCenter.mZ = (mBox[0].mZ + mBox[1].mZ) / 2;

	if (VERBOSE > 3)
		cerr << "Created residue " << mN.mResName << endl;
}

MResidue::MResidue(uint32 inNumber, char inTypeCode, MResidue* inPrevious)
	: mChainID(0)
	, mPrev(nullptr)
	, mNext(nullptr)
	, mSeqNumber(inNumber)
	, mNumber(inNumber)
	, mInsertionCode(' ')
	, mType(MapResidue(inTypeCode))
	, mSSBridgeNr(0)
	, mAccessibility(0)
	, mSecondaryStructure(loop)
	, mSheet(0)
	, mBend(false)
{
	fill(mHelixFlags, mHelixFlags + 3, helixNone);
	
	mBetaPartner[0].residue = mBetaPartner[1].residue = nullptr;
	
	mHBondDonor[0].energy = mHBondDonor[1].energy = mHBondAcceptor[0].energy = mHBondAcceptor[1].energy = 0;
	mHBondDonor[0].residue = mHBondDonor[1].residue = mHBondAcceptor[0].residue = mHBondAcceptor[1].residue = nullptr;

	static const MAtom kNullAtom = {};
	mN = mCA = mC = mO = kNullAtom;
	
	mCA.mICode = ' ';
	mCA.mResSeq = inTypeCode;
	mCA.mChainID = 'A';
}

MResidue::MResidue(const MResidue& residue)
	: mChainID(residue.mChainID)
	, mPrev(nullptr)
	, mNext(nullptr)
	, mSeqNumber(residue.mSeqNumber)
	, mNumber(residue.mNumber)
	, mType(residue.mType)
	, mSSBridgeNr(residue.mSSBridgeNr)
	, mAccessibility(residue.mAccessibility)
	, mSecondaryStructure(residue.mSecondaryStructure)
	, mC(residue.mC)
	, mN(residue.mN)
	, mCA(residue.mCA)
	, mO(residue.mO)
	, mH(residue.mH)
	, mSideChain(residue.mSideChain)
	, mSheet(residue.mSheet)
	, mBend(residue.mBend)
	, mCenter(residue.mCenter)
	, mRadius(residue.mRadius)
{
	copy(residue.mHBondDonor, residue.mHBondDonor + 2, mHBondDonor);
	copy(residue.mHBondAcceptor, residue.mHBondAcceptor + 2, mHBondAcceptor);
	copy(residue.mBetaPartner, residue.mBetaPartner + 2, mBetaPartner);
	copy(residue.mHelixFlags, residue.mHelixFlags + 3, mHelixFlags);
	copy(residue.mBox, residue.mBox + 2, mBox);
}

void MResidue::SetPrev(MResidue* inResidue)
{
	mPrev = inResidue;
	mPrev->mNext = this;
}

bool MResidue::NoChainBreak(const MResidue* from, const MResidue* to)
{
	bool result = true;
	for (const MResidue* r = from; result and r != to; r = r->mNext)
	{
		MResidue* next = r->mNext;
		if (next == nullptr)
			result = false;
		else
			result = next->mNumber == r->mNumber + 1;
	}
	return result;
}

void MResidue::SetChainID(char inID)
{
	mChainID = inID;
	
	mC.SetChainID(inID);
	mCA.SetChainID(inID);
	mO.SetChainID(inID);
	mN.SetChainID(inID);
	mH.SetChainID(inID);
	for_each(mSideChain.begin(), mSideChain.end(), boost::bind(&MAtom::SetChainID, _1, inID));
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
	if (mPrev != nullptr and NoChainBreak(mPrev, this))
		result = DihedralAngle(mPrev->GetC(), GetN(), GetCAlpha(), GetC());
	return result;
}

double MResidue::Psi() const
{
	double result = 360;
	if (mNext != nullptr and NoChainBreak(this, mNext))
		result = DihedralAngle(GetN(), GetCAlpha(), GetC(), mNext->GetN());
	return result;
}

tr1::tuple<double,char> MResidue::Alpha() const
{
	double alhpa = 360;
	char chirality = ' ';
	
	const MResidue* nextNext = mNext ? mNext->Next() : nullptr;
	if (mPrev != nullptr and nextNext != nullptr and NoChainBreak(mPrev, nextNext))
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
	const MResidue* prevPrev = mPrev ? mPrev->Prev() : nullptr;
	const MResidue* nextNext = mNext ? mNext->Next() : nullptr;
	if (prevPrev != nullptr and nextNext != nullptr and NoChainBreak(prevPrev, nextNext))
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
	if (mPrev != nullptr and NoChainBreak(mPrev, this))
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

// TODO: use the angle to improve bond energy calculation.
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

		// DSSP compatibility mode:
		result = bm::round(result * 1000) / 1000;

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
										//                                          
										//		c -> d		c     d
										
	if (a and c and NoChainBreak(a, c) and d and f and NoChainBreak(d, f))
	{
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

void MResidue::CalculateSurface(const vector<MResidue*>& inResidues)
{
	vector<MResidue*> neighbours;
	
	foreach (MResidue* r, inResidues)
	{
		MPoint center;
		double radius;
		r->GetCenterAndRadius(center, radius);
		
		if (Distance(mCenter, center) < mRadius + radius)
			neighbours.push_back(r);
	}

	mAccessibility = CalculateSurface(mN, kRadiusN, neighbours) +
					 CalculateSurface(mCA, kRadiusCA, neighbours) +
					 CalculateSurface(mC, kRadiusC, neighbours) +
					 CalculateSurface(mO, kRadiusO, neighbours);
	
	foreach (const MAtom& atom, mSideChain)
		mAccessibility += CalculateSurface(atom, kRadiusSideAtom, neighbours);
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

double MResidue::CalculateSurface(const MAtom& inAtom, double inRadius, const vector<MResidue*>& inResidues)
{
	MAccumulator accumulate;

	foreach (MResidue* r, inResidues)
	{
		if (r->AtomIntersectsBox(inAtom, inRadius))
		{
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
	
	MSurfaceDots& surfaceDots = MSurfaceDots::Instance();
	
	for (uint32 i = 0; i < surfaceDots.size(); ++i)
	{
		MPoint xx = surfaceDots[i] * radius;
		
		bool free = true;
		for (uint32 k = 0; free and k < accumulate.m_x.size(); ++k)
			free = accumulate.m_x[k].radius < DistanceSquared(xx, accumulate.m_x[k].location);
		
		if (free)
			surface += surfaceDots.weight();
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
	
	for_each(mSideChain.begin(), mSideChain.end(), boost::bind(&MAtom::WritePDB, _1, boost::ref(os)));
}

// --------------------------------------------------------------------

MChain::MChain(const MChain& chain)
	: mChainID(chain.mChainID)
{
	MResidue* previous = nullptr;
	
	foreach (const MResidue* residue, chain.mResidues)
	{
		MResidue* newResidue = new MResidue(*residue);
		newResidue->SetPrev(previous);
		mResidues.push_back(newResidue);
		previous = newResidue;
	}
}

MChain::~MChain()
{
	foreach (MResidue* residue, mResidues)
		delete residue;
}

MChain& MChain::operator=(const MChain& chain)
{
	foreach (MResidue* residue, mResidues)
		delete residue;
	mResidues.clear();

	foreach (const MResidue* residue, chain.mResidues)
		mResidues.push_back(new MResidue(*residue));
	
	mChainID = chain.mChainID;
	
	return *this;
}

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
	for_each(mResidues.begin(), mResidues.end(), boost::bind(&MResidue::WritePDB, _1, boost::ref(os)));
	
	boost::format ter("TER    %4.4d      %3.3s %c%4.4d%c");
	
	MResidue* last = mResidues.back();
	
	os << (ter % (last->GetCAlpha().mSerial + 1) % kResidueInfo[last->GetType()].name % mChainID % last->GetNumber() % ' ') << endl;
}

MResidue* MChain::GetResidueBySeqNumber(uint16 inSeqNumber, char inInsertionCode)
{
	vector<MResidue*>::iterator r = find_if(mResidues.begin(), mResidues.end(),
		boost::bind(&MResidue::GetSeqNumber, _1) == inSeqNumber and
		boost::bind(&MResidue::GetInsertionCode, _1) == inInsertionCode);
	if (r == mResidues.end())
		throw mas_exception(boost::format("Residue %d%c not found") % inSeqNumber % inInsertionCode);
	return *r;
}

void MChain::GetSequence(string& outSequence) const
{
	foreach (const MResidue* r, GetResidues())
		outSequence += kResidueInfo[r->GetType()].code;
}

// --------------------------------------------------------------------

struct MResidueID
{
	char			chain;
	uint16			seqNumber;
	char			insertionCode;
	
	bool			operator<(const MResidueID& o) const
					{
						return
							 chain < o.chain or
							(chain == o.chain and seqNumber < o.seqNumber) or
							(chain == o.chain and seqNumber == o.seqNumber and insertionCode < o.insertionCode);
					}

	bool			operator!=(const MResidueID& o) const
					{
						return chain != o.chain or seqNumber != o.seqNumber or insertionCode != o.insertionCode;
					}
};

MProtein::MProtein(istream& is, bool cAlphaOnly)
	: mResidueCount(0)
	, mChainBreaks(0)
	, mIgnoredWaterMolecules(0)
	, mNrOfHBondsInParallelBridges(0)
	, mNrOfHBondsInAntiparallelBridges(0)
{
	fill(mParallelBridgesPerLadderHistogram, mParallelBridgesPerLadderHistogram + kHistogramSize, 0);
	fill(mAntiparallelBridgesPerLadderHistogram, mAntiparallelBridgesPerLadderHistogram + kHistogramSize, 0);
	fill(mLaddersPerSheetHistogram, mLaddersPerSheetHistogram + kHistogramSize, 0);

	vector<pair<MResidueID,MResidueID> > ssbonds;
	set<char> terminatedChains;

	bool model = false;
	vector<MAtom> atoms;
	char firstAltLoc = 0;
	bool atomSeen = false;
	
	while (not is.eof())
	{
		string line;
		getline(is, line);

		if (VERBOSE > 3)
			cerr << line << endl;
		
		if (ba::starts_with(line, "HEADER"))
		{
			mHeader = line;
			ba::trim(mHeader);
			if (line.length() >= 66)
				mID = line.substr(62, 4);
			else
				mID = "UNDF";
			continue;
		}

		if (ba::starts_with(line, "COMPND"))
		{
			ba::trim(line);
			mCompound.push_back(line);
			continue;
		}

		if (ba::starts_with(line, "SOURCE"))
		{
			ba::trim(line);
			mSource.push_back(line);
			continue;
		}
		
		if (ba::starts_with(line, "AUTHOR"))
		{
			ba::trim(line);
			mAuthor.push_back(line);
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
			ssbond.first.insertionCode = line[21];
			ssbond.second.chain = line[29];
			ssbond.second.seqNumber = boost::lexical_cast<uint16>(ba::trim_copy(line.substr(30, 5)));
			ssbond.second.insertionCode = line[35];
			
			ssbonds.push_back(ssbond);
			continue;
		}
		
		// add ATOMs only if the chain isn't terminated
		if (terminatedChains.count(line[21]))
			continue;
		
		if (ba::starts_with(line, "TER   "))
		{
			if (atoms.empty())
			{
				cerr << "no atoms read before TER record " << endl
					 << line << endl;
				continue;
			}
			
			AddResidue(atoms);
			atoms.clear();
			firstAltLoc = 0;
			atomSeen = false;
			
			terminatedChains.insert(line[21]);

			continue;
		}

		if (ba::starts_with(line, "ATOM  ") or ba::starts_with(line, "HETATM"))
			//	1 - 6	Record name "ATOM "
		{
			if (cAlphaOnly and line.substr(12, 4) != " CA ")
				continue;

			atomSeen = ba::starts_with(line, "ATOM  ");

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
			
//			alternative test, check chain ID as well.
			if (not atoms.empty() and
				(atom.mChainID != atoms.back().mChainID or 
				 (atom.mResSeq != atoms.back().mResSeq or
				  (atom.mResSeq == atoms.back().mResSeq and atom.mICode != atoms.back().mICode))))
//			if (not atoms.empty() and
//				(atom.mResSeq != atoms.back().mResSeq or (atom.mResSeq == atoms.back().mResSeq and atom.mICode != atoms.back().mICode)))
			{
				AddResidue(atoms);
				atoms.clear();
				firstAltLoc = 0;
			}

			try
			{
				atom.mType = MapElement(line.substr(76, 2));
			}
			catch (exception& e)
			{
				if (VERBOSE)
					cerr << e.what() << endl;
				atom.mType = kUnknownAtom;
			}

			if (atom.mType == kHydrogen)
				continue;
			
			if (atom.mAltLoc != ' ')
			{
				if (firstAltLoc == 0)
					firstAltLoc = atom.mAltLoc;
				if (atom.mAltLoc == firstAltLoc)
					atom.mAltLoc = 'A';
			}

			if (firstAltLoc != 0 and atom.mAltLoc != ' ' and atom.mAltLoc != firstAltLoc)
			{
				if (VERBOSE)
					cerr << "skipping alternate atom record " << atom.mResName << endl;
				continue;
			}
			
			atoms.push_back(atom);
		}
	}
	
	if (not atoms.empty())	// we have read atoms without a TER
	{
		if (atomSeen and VERBOSE)
			cerr << "ATOM records not terminated by TER record" << endl;
		
		AddResidue(atoms);
	}

	// map the sulfur bridges
	uint32 ssbondNr = 1;
	typedef pair<MResidueID,MResidueID> SSBond;
	foreach (const SSBond& ssbond, ssbonds)
	{
		try
		{
			MResidue* first = GetResidue(ssbond.first.chain, ssbond.first.seqNumber, ssbond.first.insertionCode);
			MResidue* second = GetResidue(ssbond.second.chain, ssbond.second.seqNumber, ssbond.second.insertionCode);
		
			if (first == second)
				throw mas_exception("first and second residue are the same");
		
			first->SetSSBridgeNr(ssbondNr);
			second->SetSSBridgeNr(ssbondNr);
			
			mSSBonds.push_back(make_pair(first, second));
			
			++ssbondNr;
		}
		catch (exception& e)
		{
			if (VERBOSE)
				cerr << "invalid residue referenced in SSBOND record: " << e.what() << endl;
		}
	}
	
	mChains.erase(
		remove_if(mChains.begin(), mChains.end(), boost::bind(&MChain::Empty, _1)),
		mChains.end());

	if (VERBOSE and mIgnoredWaterMolecules)
		cerr << "Ignored " << mIgnoredWaterMolecules << " water molecules" << endl;
	
	if (mChains.empty())
		throw mas_exception("empty protein, or no valid complete residues");
}

MProtein::MProtein(const string& inID, MChain* inChain)
	: mID(inID)
{
	mChains.push_back(inChain);
}

MProtein::~MProtein()
{
	foreach (MChain* chain, mChains)
		delete chain;
}

string MProtein::GetCompound() const
{
	string result;
	if (not mCompound.empty())
	{
		result = mCompound.front();
		if (ba::starts_with(result.substr(10), "MOL_ID: "))
		{
			if (mCompound.size() > 1)
				result = mCompound[1];
			else
				result.clear();
		}
	}
	return result;
}

string MProtein::GetSource() const
{
	string result;
	if (not mSource.empty())
	{
		result = mSource.front();
		if (ba::starts_with(result.substr(10), "MOL_ID: "))
		{
			if (mSource.size() > 1)
				result = mSource[1];
			else
				result.clear();
		}
	}
	return result;
}

string MProtein::GetAuthor() const
{
	string result;
	if (not mAuthor.empty())
	{
		result = mAuthor.front();
		if (ba::starts_with(result.substr(10), "MOL_ID: "))
		{
			if (mAuthor.size() > 1)
				result = mAuthor[1];
			else
				result.clear();
		}
	}
	return result;
}

void MProtein::GetStatistics(uint32& outNrOfResidues, uint32& outNrOfChains,
	uint32& outNrOfSSBridges, uint32& outNrOfIntraChainSSBridges,
	uint32& outNrOfHBonds, uint32 outNrOfHBondsPerDistance[11]) const
{
	outNrOfResidues = mResidueCount;
	outNrOfChains = mChains.size() + mChainBreaks;
	outNrOfSSBridges = mSSBonds.size();
	
	outNrOfIntraChainSSBridges = 0;
	for (vector<pair<MResidue*,MResidue*> >::const_iterator ri = mSSBonds.begin(); ri != mSSBonds.end(); ++ri)
	{
		if (ri->first->GetChainID() == ri->second->GetChainID() and
			(MResidue::NoChainBreak(ri->first, ri->second) or MResidue::NoChainBreak(ri->first, ri->second)))
		{
			++outNrOfIntraChainSSBridges;
		}
	}
	
	outNrOfHBonds = 0;
	foreach (const MChain* chain, mChains)
	{
		foreach (const MResidue* r, chain->GetResidues())
		{
			const HBond* donor = r->Donor();
			
			for (uint32 i = 0; i < 2; ++i)
			{
				if (donor[i].residue != nullptr and donor[i].energy < kMaxHBondEnergy)
				{
					++outNrOfHBonds;
					int32 k = donor[i].residue->GetNumber() - r->GetNumber();
					if (k >= -5 and k <= 5)
						outNrOfHBondsPerDistance[k + 5] += 1;
				}
			}
		}
	}
}

void MProtein::GetResiduesPerAlphaHelixHistogram(uint32 outHistogram[30]) const
{
	fill(outHistogram, outHistogram + 30, 0);

	foreach (const MChain* chain, mChains)
	{
		uint32 helixLength = 0;
		
		foreach (const MResidue* r, chain->GetResidues())
		{
			if (r->GetSecondaryStructure() == alphahelix)
				++helixLength;
			else if (helixLength > 0)
			{
				if (helixLength > kHistogramSize)
					helixLength = kHistogramSize;
				
				outHistogram[helixLength - 1] += 1;
				helixLength = 0;
			}
		}
	}
}

void MProtein::GetParallelBridgesPerLadderHistogram(uint32 outHistogram[30]) const
{
	copy(mParallelBridgesPerLadderHistogram, mParallelBridgesPerLadderHistogram + kHistogramSize, outHistogram);
}

void MProtein::GetAntiparallelBridgesPerLadderHistogram(uint32 outHistogram[30]) const
{
	copy(mAntiparallelBridgesPerLadderHistogram, mAntiparallelBridgesPerLadderHistogram + kHistogramSize, outHistogram);
}

void MProtein::GetLaddersPerSheetHistogram(uint32 outHistogram[30]) const
{
	copy(mLaddersPerSheetHistogram, mLaddersPerSheetHistogram + kHistogramSize, outHistogram);
}

void MProtein::AddResidue(const vector<MAtom>& inAtoms)
{
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
		MChain& chain = GetChain(inAtoms.front().mChainID);
		vector<MResidue*>& residues(chain.GetResidues());
	
		MResidue* prev = nullptr;
		if (not residues.empty())
			prev = residues.back();

		uint32 resNumber = mResidueCount + mChains.size() + mChainBreaks;
		MResidue* r = new MResidue(resNumber, prev, inAtoms);
		if (prev != nullptr and not prev->ValidDistance(*r))	// check for chain breaks
		{
			if (VERBOSE)
				cerr << boost::format("The distance between residue %1% and %2% is larger than the maximum peptide bond length")
						% prev->GetNumber() % resNumber << endl;
			
			++mChainBreaks;
			r->SetNumber(resNumber + 1);
		}
		
		residues.push_back(r);
		++mResidueCount;
	}
	else if (string(inAtoms.front().mResName) == "HOH ")
		++mIgnoredWaterMolecules;
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
	vector<MResidue*> residues;
	residues.reserve(mResidueCount);
	foreach (const MChain* chain, mChains)
		residues.insert(residues.end(), chain->GetResidues().begin(), chain->GetResidues().end());
	
	if (VERBOSE)
		cerr << "using " << residues.size() << " residues" << endl;

	boost::thread t(boost::bind(&MProtein::CalculateAccessibilities, this, boost::ref(residues)));

	CalculateHBondEnergies(residues);
	CalculateBetaSheets(residues);
	CalculateAlphaHelices(residues);

	t.join();
}

void MProtein::CalculateHBondEnergies(const std::vector<MResidue*>& inResidues)
{
	if (VERBOSE)
		cerr << "Calculate H-bond energies" << endl;
	
	// Calculate the HBond energies
	for (uint32 i = 0; i + 1 < inResidues.size(); ++i)
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

// TODO: improve alpha helix calculation by better recognizing pi-helices 
void MProtein::CalculateAlphaHelices(const std::vector<MResidue*>& inResidues)
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
			
			for (uint32 i = 0; i + stride < res.size(); ++i)
			{
				if (MResidue::TestBond(res[i + stride], res[i]) and MResidue::NoChainBreak(res[i], res[i + stride]))
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

	for (uint32 i = 1; i + 4 < inResidues.size(); ++i)
	{
		if (inResidues[i]->IsHelixStart(4) and inResidues[i - 1]->IsHelixStart(4))
		{
			for (uint32 j = i; j <= i + 3; ++j)
				inResidues[j]->SetSecondaryStructure(alphahelix);
		}
	}

	for (uint32 i = 1; i + 3 < inResidues.size(); ++i)
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

	for (uint32 i = 1; i + 5 < inResidues.size(); ++i)
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
			
	for (uint32 i = 1; i + 1 < inResidues.size(); ++i)
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

void MProtein::CalculateBetaSheets(const std::vector<MResidue*>& inResidues)
{
	if (VERBOSE)
		cerr << "Calculate beta sheets" << endl;
	
	// Calculate Bridges
	vector<MBridge> bridges;
	if (inResidues.size() > 4)
	{
		for (uint32 i = 1; i + 4 < inResidues.size(); ++i)
		{
			MResidue* ri = inResidues[i];
			
			for (uint32 j = i + 3; j + 1 < inResidues.size(); ++j)
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
			uint32 ibi = bridges[i].i.front();
			uint32 iei = bridges[i].i.back();
			uint32 jbi = bridges[i].j.front();
			uint32 jei = bridges[i].j.back();
			uint32 ibj = bridges[j].i.front();
			uint32 iej = bridges[j].i.back();
			uint32 jbj = bridges[j].j.front();
			uint32 jej = bridges[j].j.back();

			if (bridges[i].type != bridges[j].type or
				MResidue::NoChainBreak(inResidues[min(ibi, ibj)], inResidues[max(iei, iej)]) == false or
				MResidue::NoChainBreak(inResidues[min(jbi, jbj)], inResidues[max(jei, jej)]) == false or
				ibj - iei >= 6 or
				(iei >= ibj and ibi <= iej))
			{
				continue;
			}
			
			bool bulge;
			if (bridges[i].type == btParallel)
				bulge = ((jbj - jei < 6 and ibj - iei < 3) or (jbj - jei < 3));
			else
				bulge = ((jbi - jej < 6 and ibj - iei < 3) or (jbi - jej < 3));

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
	{
		ladderset.insert(&bridge);
		
		uint32 n = bridge.i.size();
		if (n > kHistogramSize)
			n = kHistogramSize;
		
		if (bridge.type == btParallel)
			mParallelBridgesPerLadderHistogram[n - 1] += 1;
		else
			mAntiparallelBridgesPerLadderHistogram[n - 1] += 1;
	}
	
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
		
		uint32 nrOfLaddersPerSheet = sheetset.size();
		if (nrOfLaddersPerSheet > kHistogramSize)
			nrOfLaddersPerSheet = kHistogramSize;
		if (nrOfLaddersPerSheet == 1 and (*sheetset.begin())->i.size() > 1)
			mLaddersPerSheetHistogram[0] += 1;
		else if (nrOfLaddersPerSheet > 1)
			mLaddersPerSheetHistogram[nrOfLaddersPerSheet - 1] += 1;
		
		++sheet;
	}

	foreach (MBridge& bridge, bridges)
	{
		// find out if any of the i and j set members already have
		// a bridge assigned, if so, we're assigning bridge 2
		
		uint32 betai = 0, betaj = 0;
		
		foreach (uint32 l, bridge.i)
		{
			if (inResidues[l]->GetBetaPartner(0).residue != nullptr)
			{
				betai = 1;
				break;
			}
		}

		foreach (uint32 l, bridge.j)
		{
			if (inResidues[l]->GetBetaPartner(0).residue != nullptr)
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
			mNrOfHBondsInParallelBridges += bridge.i.back() - bridge.i.front() + 2;
			
			deque<uint32>::iterator j = bridge.j.begin();
			foreach (uint32 i, bridge.i)
				inResidues[i]->SetBetaPartner(betai, inResidues[*j++], bridge.ladder, true);

			j = bridge.i.begin();
			foreach (uint32 i, bridge.j)
				inResidues[i]->SetBetaPartner(betaj, inResidues[*j++], bridge.ladder, true);
		}
		else
		{
			mNrOfHBondsInAntiparallelBridges += bridge.i.back() - bridge.i.front() + 2;

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

void MProtein::CalculateAccessibilities(const std::vector<MResidue*>& inResidues)
{
	if (VERBOSE)
		cerr << "Calculate accessibilities" << endl;

	uint32 nr_of_threads = boost::thread::hardware_concurrency();
	if (nr_of_threads <= 1)
	{
		foreach (MResidue* residue, inResidues)
			residue->CalculateSurface(inResidues);
	}
	else
	{
		MResidueQueue queue;
	
		boost::thread_group t;
	
		for (uint32 ti = 0; ti < nr_of_threads; ++ti)
			t.create_thread(boost::bind(&MProtein::CalculateAccessibility, this,
				boost::ref(queue), boost::ref(inResidues)));
		
		foreach (MResidue* residue, inResidues)
			queue.put(residue);
		
		queue.put(nullptr);
		
		t.join_all();
	}
}

void MProtein::CalculateAccessibility(MResidueQueue& inQueue,
	const std::vector<MResidue*>& inResidues)
{
	// make sure the MSurfaceDots is constructed once
	(void)MSurfaceDots::Instance();
	
	for (;;)
	{
		MResidue* residue = inQueue.get();
		if (residue == nullptr)
			break;
		
		residue->CalculateSurface(inResidues);
	}
	
	inQueue.put(nullptr);
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

MResidue* MProtein::GetResidue(char inChainID, uint16 inSeqNumber, char inInsertionCode)
{
	MChain& chain = GetChain(inChainID);
	if (chain.GetResidues().empty())
		throw mas_exception(boost::format("Invalid chain id '%c'") % inChainID);
	return chain.GetResidueBySeqNumber(inSeqNumber, inInsertionCode);
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

void MProtein::GetSequence(char inChain, entry& outEntry) const
{
	if (inChain == 0)
		inChain = mChains.front()->GetChainID();
	
	string seq;
	foreach (const MResidue* r, GetChain(inChain).GetResidues())
	{
		seq += kResidueInfo[r->GetType()].code;
		outEntry.m_positions.push_back(r->GetSeqNumber());
		outEntry.m_ss += r->GetSecondaryStructure();
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
