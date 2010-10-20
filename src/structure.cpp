// structure related stuff

#include "MRS.h"
#include "mas.h"

#include <set>

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
	, mSecondaryStructure(loop)
	, mSheet(0)
{
	if (mPrev != nil)
		mPrev->mNext = this;
	
	mHBondDonor[0].energy = mHBondDonor[1].energy = mHBondAcceptor[0].energy = mHBondAcceptor[1].energy = 0;
	mHBondDonor[0].residue = mHBondDonor[1].residue = mHBondAcceptor[0].residue = mHBondAcceptor[1].residue = nil;
	mBetaPartner[0] = mBetaPartner[1] = nil;
	
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
	const MResidue* a = Prev();			//		  \			  /
	const MResidue* b = this;			//		b	e		b	e
	const MResidue* c = Next();			// 		  /			  \                      ..
	const MResidue* d = test->Prev();	//		c	f		c	f
	const MResidue* e = test;			//
	const MResidue* f = test->Next();	// III.	a <- f	IV. a	  f		antiparallel
										//		                                   
	MBridgeType result = nobridge;		//		b	 e      b <-> e                  
	if (a and c and d and f)			//                                          
	{									//		c -> d		c     d
		if ((TestBond(c, e) and TestBond(e, a)) or (TestBond(f, b) and TestBond(b, d)))
			result = parallel;
		else if ((TestBond(c, d) and TestBond(f, a)) or (TestBond(e, b) and TestBond(b, e)))
			result = antiparallel;
	}
	
	return result;
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

void MResidue::WriteDSSP(ostream& os)
{
/*   
  #  RESIDUE AA STRUCTURE BP1 BP2  ACC     N-H-->O    O-->H-N    N-H-->O    O-->H-N    TCO  KAPPA ALPHA  PHI   PSI    X-CA   Y-CA   Z-CA 
 */
	boost::format kDSSPResidueLine(
	"%5.5d%5.5d %c %c  %c     %c  %4.4d%4.4d%c     %11s%11s%11s%11s  %6.3f%6.1f%6.1f%6.1f%6.1f %6.1f %6.1f %6.1f");
	
	const MAtom& ca = GetCAlpha();
	
	char code = kResidueInfo[mType].code;
	if (mType == kCysteine and GetSSBridgeNr() != 0)
		code = 'a' - 1 + (GetSSBridgeNr() % 26);

	double alpha;
	char chirality;
	tr1::tie(alpha,chirality) = Alpha();
	
	char ss;
	switch (mSecondaryStructure)
	{
		case alphahelix:	ss = 'H'; break;
		case betabridge:	ss = 'B'; break;
		case strand:		ss = 'E'; break;
		case helix_3:		ss = 'G'; break;
		case helix_5:		ss = 'I'; break;
		case turn:			ss = 'T'; break;
		case bend:			ss = 'S'; break;
		case loop:			ss = ' '; break;
	}
	
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
	
	uint32 bp1 = 0, bp2 = 0;
	if (mBetaPartner[0] != nil)
		bp1 = mBetaPartner[0]->mNumber;
	if (mBetaPartner[1] != nil)
		bp2 = mBetaPartner[1]->mNumber;
	
	char sheet = ' ';
	if (mSheet != 0)
		sheet = 'A' + (mSheet - 1) % 26;
	
	cout << (kDSSPResidueLine % mNumber % ca.mResSeq % ca.mChainID % code % ss % chirality % 
		bp1 % bp2 % sheet %
		NHO1 % ONH1 % NHO2 % ONH2 %
		TCO() % Kappa() % alpha % Phi() % Psi() % ca.mLoc.mX % ca.mLoc.mY % ca.mLoc.mZ) << endl;
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

void MChain::WriteDSSP(std::ostream& os)
{
	for_each(mResidues.begin(), mResidues.end(), boost::bind(&MResidue::WriteDSSP, _1, ref(os)));
}

MResidue& MChain::GetResidueBySeqNumber(uint16 inSeqNumber)
{
	vector<MResidue*>::iterator r = find_if(mResidues.begin(), mResidues.end(),
		boost::bind(&MResidue::GetSeqNumber, _1) == inSeqNumber);
	if (r == mResidues.end())
		throw mas_exception(boost::format("Residue %d not found") % inSeqNumber);
	return **r;
}

void MChain::AddResidue(uint32 inNumber, const vector<MAtom>& inAtoms)
{
	MResidue* prev = nil;
	if (not mResidues.empty())
		prev = mResidues.back();
	mResidues.push_back(new MResidue(*this, inNumber, prev, inAtoms));
}

// --------------------------------------------------------------------

MProtein::MProtein(istream& is, bool cAlphaOnly)
{
	bool model = false;
	uint32 resNumber = 1;
	vector<MAtom> atoms;
	
	while (not is.eof())
	{
		string line;
		getline(is, line);
		
		if (ba::starts_with(line, "HEADER"))
		{
			mID = line.substr(62, 4);
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
		
		if (ba::starts_with(line, "ATOM  "))
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
			
			atom.mType = MapElement(line.substr(77, 2));
			
			if (not atoms.empty() and atom.mResSeq != atoms.back().mResSeq)
			{
				MChain& chain = GetChain(atoms.back().mChainID);
				
				chain.AddResidue(resNumber, atoms);
				++resNumber;
				
				atoms.clear();
			}
			
			atoms.push_back(atom);
		}
	}

	if (not atoms.empty())
	{
		MChain& chain = GetChain(atoms.back().mChainID);
		chain.AddResidue(resNumber, atoms);

		foreach (const MResidue* r, GetChain(atoms.back().mChainID).GetResidues())
		{
			cerr << r->GetNumber() << '\t' << r->GetSeqNumber() << '\t' << kResidueInfo[r->GetType()].code << endl;
		}
	}
}

MChain& MProtein::GetChain(char inChainID)
{
	map<char,MChain>::iterator chain = mChains.find(inChainID);
	if (chain == mChains.end())
		mChains[inChainID] = MChain(inChainID);
	return mChains[inChainID];
}

void MProtein::GetPoints(std::vector<MPoint>& outPoints) const
{
	for (std::map<char,MChain>::const_iterator chain = mChains.begin(); chain != mChains.end(); ++chain)
	{
		foreach (const MResidue* r, chain->second.GetResidues())
			r->GetPoints(outPoints);
	}
}

void MProtein::Translate(const MPoint& inTranslation)
{
	for (map<char,MChain>::iterator chain = mChains.begin(); chain != mChains.end(); ++chain)
		chain->second.Translate(inTranslation);
}

void MProtein::Rotate(const MQuaternion& inRotation)
{
	for (map<char,MChain>::iterator chain = mChains.begin(); chain != mChains.end(); ++chain)
		chain->second.Rotate(inRotation);
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
	
	if (VERBOSE)
		cerr << "Calculate HBond energies" << endl;
	
	// Calculate the HBond energies
	vector<MResidue*> residues;
	for (map<char,MChain>::iterator chain = mChains.begin(); chain != mChains.end(); ++chain)
		residues.insert(residues.end(), chain->second.GetResidues().begin(), chain->second.GetResidues().end());
	
	for (uint32 i = 0; i < residues.size() - 1; ++i)
	{
		MResidue* ri = residues[i];
		
		for (uint32 j = i + 1; j < residues.size(); ++j)
		{
			MResidue* rj = residues[j];
			
			if (Distance(ri->GetCAlpha(), rj->GetCAlpha()) < kMinimalCADistance)
			{
				MResidue::CalculateHBondEnergy(*ri, *rj);
				if (j != i + 1)
					MResidue::CalculateHBondEnergy(*rj, *ri);
			}
		}
	}

	if (VERBOSE)
		cerr << "Calculate Bridges" << endl;
	
	// Calculate Bridges
	vector<MBridge> bridges;
	for (uint32 i = 1; i < residues.size() - 4; ++i)
	{
		MResidue* ri = residues[i];
		
		for (uint32 j = i + 3; j < residues.size() - 1; ++j)
		{
			MResidue* rj = residues[j];
			
			MBridgeType type = ri->TestBridge(rj);
			if (type == nobridge)
				continue;
			
			bool found = false;
			foreach (MBridge& bridge, bridges)
			{
				if (type != bridge.type or i != bridge.i.back() + 1)
					continue;
				
				if (type == parallel and bridge.j.back() + 1 == j)
				{
					bridge.i.push_back(i);
					bridge.j.push_back(j);
					found = true;
					break;
				}

				if (type == antiparallel and bridge.j.front() - 1 == j)
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
			if (bridges[i].type == parallel)
				bulge = (jbj - jei < 6 and ibj - iei < 3) or (jbj - jei < 3);
			else
				bulge = (jbi - jej < 6 and ibj - iei < 3) or (jbi - jej < 3);

			if (bulge)
			{
				bridges[i].i.insert(bridges[i].i.end(), bridges[j].i.begin(), bridges[j].i.end());
				if (bridges[i].type == parallel)
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
	
	uint32 sheet = 1, ladder = 1;
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
			if (residues[l]->GetBetaPartner(0) != nil)
			{
				betai = 1;
				break;
			}
		}

		foreach (uint32 l, bridge.j)
		{
			if (residues[l]->GetBetaPartner(0) != nil)
			{
				betaj = 1;
				break;
			}
		}
		
		MSecondaryStructure ss = betabridge;
		if (bridge.i.size() > 1)
			ss = strand;
		
		if (bridge.type == parallel)
		{
			deque<uint32>::iterator j = bridge.j.begin();

			foreach (uint32 i, bridge.i)
			{
				residues[i]->SetBetaPartner(betai, residues[*j++], bridge.ladder);
				residues[i]->SetSheet(bridge.sheet);
			}

			j = bridge.i.begin();

			foreach (uint32 i, bridge.j)
			{
				residues[i]->SetBetaPartner(betaj, residues[*j++], bridge.ladder);
				residues[i]->SetSheet(bridge.sheet);
			}
		}
		else
		{
			deque<uint32>::reverse_iterator j = bridge.j.rbegin();

			foreach (uint32 i, bridge.i)
			{
				residues[i]->SetBetaPartner(betai, residues[*j++], bridge.ladder);
				residues[i]->SetSheet(bridge.sheet);
			}

			j = bridge.i.rbegin();

			foreach (uint32 i, bridge.j)
			{
				residues[i]->SetBetaPartner(betaj, residues[*j++], bridge.ladder);
				residues[i]->SetSheet(bridge.sheet);
			}
		}

		for (uint32 i = bridge.i.front(); i <= bridge.i.back(); ++i)
		{
			if (residues[i]->GetSecondaryStructure() != strand)
				residues[i]->SetSecondaryStructure(ss);
			residues[i]->SetSheet(bridge.sheet);
		}

		for (uint32 i = bridge.j.front(); i <= bridge.j.back(); ++i)
		{
			if (residues[i]->GetSecondaryStructure() != strand)
				residues[i]->SetSecondaryStructure(ss);
			residues[i]->SetSheet(bridge.sheet);
		}
	}
}

void MProtein::Center()
{
	vector<MPoint> p;
	GetPoints(p);
	
	MPoint t = CenterPoints(p);
	
	Translate(MPoint(-t.mX, -t.mY, -t.mZ));
}

MResidue& MProtein::GetResidue(char inChainID, uint16 inSeqNumber)
{
	MChain& chain = mChains[inChainID];
	if (chain.GetResidues().empty())
		throw mas_exception(boost::format("Invalid chain id '%c'") % inChainID);
	return chain.GetResidueBySeqNumber(inSeqNumber);
}

void MProtein::GetCAlphaLocations(char inChain, vector<MPoint>& outPoints) const
{
	if (inChain == 0)
		inChain = mChains.begin()->first;
	
	map<char,MChain>::const_iterator chain = mChains.find(inChain);
	assert(chain != mChains.end());
	
	foreach (const MResidue* r, chain->second.GetResidues())
		outPoints.push_back(r->GetCAlpha());
}

MPoint MProtein::GetCAlphaPosition(char inChain, int16 inPDBResSeq) const
{
	if (inChain == 0)
		inChain = mChains.begin()->first;
	
	map<char,MChain>::const_iterator chain = mChains.find(inChain);
	assert(chain != mChains.end());
	
	MPoint result;
	
	foreach (const MResidue* r, chain->second.GetResidues())
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
		inChain = mChains.begin()->first;
	
	map<char,MChain>::const_iterator chain = mChains.find(inChain);
	assert(chain != mChains.end());
	
	string seq;
	
	foreach (const MResidue* r, chain->second.GetResidues())
	{
		seq += kResidueInfo[r->GetType()].code;
		outEntry.m_positions.push_back(r->GetSeqNumber());
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
	
	foreach (const MResidue* r, chain->second.GetResidues())
		seq += kResidueInfo[r->GetType()].code;
	
	outSequence = encode(seq);
}

void MProtein::WritePDB(ostream& os)
{
	for (map<char,MChain>::iterator c = mChains.begin(); c != mChains.end(); ++c)
		c->second.WritePDB(os);
}

void MProtein::WriteDSSP(ostream& os)
{
	os << "  #  RESIDUE AA STRUCTURE BP1 BP2  ACC     N-H-->O    O-->H-N    N-H-->O    O-->H-N    TCO  KAPPA ALPHA  PHI   PSI    X-CA   Y-CA   Z-CA " << endl;
	boost::format kDSSPResidueLine(
	"%5.5d        !*             0   0    0      0, 0.0     0, 0.0     0, 0.0     0, 0.0   0.000 360.0 360.0 360.0 360.0    0.0    0.0    0.0");  

	for (map<char,MChain>::iterator c = mChains.begin(); c != mChains.end(); ++c)
	{
		c->second.WriteDSSP(os);
		if (next(c) != mChains.end())
			os << (kDSSPResidueLine % (c->second.GetResidues().back()->GetNumber() + 1)) << endl;
	}
}
