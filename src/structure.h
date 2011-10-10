// Copyright Maarten L. Hekkelman, Radboud University 2008-2011.
//   Distributed under the Boost Software License, Version 1.0.
//       (See accompanying file LICENSE_1_0.txt or copy at    
//             http://www.boost.org/LICENSE_1_0.txt)      

#pragma once

#include "primitives-3d.h"
#include "align-2d.h"

struct MAtom;
class MResidue;
class MChain;
class MProtein;

// forward declaration of buffer
template<typename T, uint32 N> class buffer;
typedef buffer<MResidue*,100>	MResidueQueue;

const uint32 kHistogramSize = 30;

// a limited set of known atoms. This is an obvious candidate for improvement of DSSP.
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
	kMagnesium,
	kPotassium,
	kCalcium,
	kZinc,
	kSelenium,
	
	kAtomTypeCount
};

MAtomType MapElement(std::string inElement);

// for now, MAtom contains exactly what the ATOM line contains in a PDB file
struct MAtom
{
	uint32		mSerial;
	char		mName[5];
	char		mAltLoc;
	char		mResName[5];
	char		mChainID;
	int16		mResSeq;
	char		mICode;
	MAtomType	mType;
	MPoint		mLoc;
	double		mOccupancy;
	double		mTempFactor;
	char		mElement[3];
	int			mCharge;

	void		SetChainID(char inID)					{ mChainID = inID;}
	std::string	GetName() const							{ return mName; }
	void		Translate(const MPoint& inTranslation)	{ mLoc += inTranslation; }
	void		Rotate(const MQuaternion& inRotation)	{ mLoc.Rotate(inRotation); }
	void		WritePDB(std::ostream& os) const;

				operator const MPoint&() const			{ return mLoc; }
				operator MPoint&()						{ return mLoc; }
};

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
};

// a residue number to info mapping
extern const MResidueInfo kResidueInfo[];

MResidueType MapResidue(std::string inName);

struct HBond
{
	MResidue*		residue;
	double			energy;
};

enum MBridgeType
{
	btNoBridge, btParallel, btAntiParallel
};

struct MBridgeParner
{
	MResidue*		residue;
	uint32			ladder;
	bool			parallel;
};

enum MHelixFlag
{
	helixNone, helixStart, helixEnd, helixStartAndEnd, helixMiddle
};

enum MSecondaryStructure
{
	loop,		//' '
	alphahelix,	// H
	betabridge, // B
	strand,		// E
	helix_3,	// G
	helix_5,	// I
	turn,		// T
	bend		// S
};

class MResidue
{
  public:
						MResidue(const MResidue& residue);
						MResidue(uint32 inNumber, char inTypeCode, MResidue* inPrevious);
						MResidue(uint32 inNumber,
							MResidue* inPrevious, const std::vector<MAtom>& inAtoms);

	void				SetChainID(char inID);
	char				GetChainID() const				{ return mChainID; }

	MResidueType		GetType() const					{ return mType; }

	const MAtom&		GetCAlpha() const				{ return mCA; }
	const MAtom&		GetC() const					{ return mC; }
	const MAtom&		GetN() const					{ return mN; }
	const MAtom&		GetO() const					{ return mO; }
	const MAtom&		GetH() const					{ return mH; }

	double				Phi() const;
	double				Psi() const;
	std::tr1::tuple<double,char>
						Alpha() const;
	double				Kappa() const;
	double				TCO() const;
	
	double				Accessibility() const			{ return mAccessibility; }
	
	void				SetSecondaryStructure(MSecondaryStructure inSS)
														{ mSecondaryStructure = inSS; }
	MSecondaryStructure	GetSecondaryStructure() const	{ return mSecondaryStructure; }
	
	const MResidue*		Next() const					{ return mNext; }
	const MResidue*		Prev() const					{ return mPrev; }

	void				SetPrev(MResidue* inResidue);
	
	void				SetBetaPartner(uint32 n, MResidue* inResidue, uint32 inLadder,
							bool inParallel);
	MBridgeParner		GetBetaPartner(uint32 n) const;
						
	void				SetSheet(uint32 inSheet)	{ mSheet = inSheet; }
	uint32				GetSheet() const			{ return mSheet; }
	
	bool				IsBend() const				{ return mBend; }
	void				SetBend(bool inBend)		{ mBend = inBend; }
	
	MHelixFlag			GetHelixFlag(uint32 inHelixStride) const;
	bool				IsHelixStart(uint32 inHelixStride) const;
	void				SetHelixFlag(uint32 inHelixStride, MHelixFlag inHelixFlag);

	void				SetSSBridgeNr(uint8 inBridgeNr);
	uint8				GetSSBridgeNr() const;

	void				AddAtom(MAtom& inAtom);
	
	HBond*				Donor()						{ return mHBondDonor; }
	HBond*				Acceptor()					{ return mHBondAcceptor; }

	const HBond*		Donor() const				{ return mHBondDonor; }
	const HBond*		Acceptor() const			{ return mHBondAcceptor; }

	bool				ValidDistance(const MResidue& inNext) const;

	static bool			TestBond(const MResidue* a, const MResidue* b)
						{
							return a->TestBond(b);
						}

	// bridge functions
	MBridgeType			TestBridge(MResidue* inResidue) const;

	uint16				GetSeqNumber() const		{ return mSeqNumber; }
	char				GetInsertionCode() const	{ return mInsertionCode; }
	
	void				SetNumber(uint16 inNumber)	{ mNumber = inNumber; }
	uint16				GetNumber() const			{ return mNumber; }

	void				Translate(const MPoint& inTranslation);
	void				Rotate(const MQuaternion& inRotation);

	void				WritePDB(std::ostream& os);

	static double		CalculateHBondEnergy(MResidue& inDonor, MResidue& inAcceptor);

	std::vector<MAtom>&	GetSideChain()				{ return mSideChain; }
	const std::vector<MAtom>&
						GetSideChain() const		{ return mSideChain; }

	void				GetPoints(std::vector<MPoint>& outPoints) const;

	void				CalculateSurface(const std::vector<MResidue*>& inResidues);

	void				GetCenterAndRadius(MPoint& outCenter, double& outRadius) const
													{ outCenter = mCenter; outRadius = mRadius; }

	static bool			NoChainBreak(const MResidue* from, const MResidue* to);

  protected:

	double				CalculateSurface(
							const MAtom& inAtom, double inRadius,
							const std::vector<MResidue*>& inResidues);

	bool				TestBond(const MResidue* other) const;

	void				ExtendBox(const MAtom& atom, double inRadius);
	bool				AtomIntersectsBox(const MAtom& atom, double inRadius) const;

	char				mChainID;
	MResidue*			mPrev;
	MResidue*			mNext;
	int32				mSeqNumber, mNumber;
	char				mInsertionCode;
	MResidueType		mType;
	uint8				mSSBridgeNr;
	double				mAccessibility;
	MSecondaryStructure	mSecondaryStructure;
	MAtom				mC, mN, mCA, mO, mH;
	HBond				mHBondDonor[2], mHBondAcceptor[2];
	std::vector<MAtom>	mSideChain;
	MBridgeParner		mBetaPartner[2];
	uint32				mSheet;
	MHelixFlag			mHelixFlags[3];	//
	bool				mBend;
	MPoint				mBox[2];		// The 3D box containing all atoms
	MPoint				mCenter;		// and the 3d Sphere containing all atoms
	double				mRadius;

  private:
	MResidue&			operator=(const MResidue& residue);
};

class MChain
{
  public:

						MChain(const MChain& chain);
						MChain(char inChainID = 0) : mChainID(inChainID) {}
						~MChain();

	MChain&				operator=(const MChain& chain);

	char				GetChainID() const					{ return mChainID; }
	void				SetChainID(char inID);

	MResidue*			GetResidueBySeqNumber(uint16 inSeqNumber, char inInsertionCode);
	
	void				GetSequence(std::string& outSequence) const;

	void				Translate(const MPoint& inTranslation);
	void				Rotate(const MQuaternion& inRotation);

	void				WritePDB(std::ostream& os);
	
	std::vector<MResidue*>&
						GetResidues()						{ return mResidues; }
	const std::vector<MResidue*>&
						GetResidues() const					{ return mResidues; }

	bool				Empty() const						{ return mResidues.empty(); }

  private:
	char				mChainID;
	std::vector<MResidue*>
						mResidues;
};

class MProtein
{
  public:
						MProtein() {}
						MProtein(const std::string& inID, MChain* inChain);
						~MProtein();

	const std::string&	GetID() const					{ return mID; }
						
						MProtein(std::istream& is, bool inCAlphaOnly = false);
	
	const std::string&	GetHeader() const				{ return mHeader; }
	std::string			GetCompound() const;
	std::string			GetSource() const;
	std::string			GetAuthor() const;

	void				CalculateSecondaryStructure();
	
	void				GetStatistics(uint32& outNrOfResidues, uint32& outNrOfChains,
							uint32& outNrOfSSBridges, uint32& outNrOfIntraChainSSBridges,
							uint32& outNrOfHBonds, uint32 outNrOfHBondsPerDistance[11]) const;
	
	void				GetCAlphaLocations(char inChain, std::vector<MPoint>& outPoints) const;
	MPoint				GetCAlphaPosition(char inChain, int16 inPDBResSeq) const;
	
	void				GetSequence(char inChain, entry& outEntry) const;
	void				GetSequence(char inChain, sequence& outSequence) const;

	void				Center();
	void				Translate(const MPoint& inTranslation);
	void				Rotate(const MQuaternion& inRotation);

	void				WritePDB(std::ostream& os);
	
	void				GetPoints(std::vector<MPoint>& outPoints) const;

	char				GetFirstChainID() const								{ return mChains.front()->GetChainID(); }

	void				SetChain(char inChainID, const MChain& inChain);

	MChain&				GetChain(char inChainID);
	const MChain&		GetChain(char inChainID) const;
	
	const std::vector<MChain*>&
						GetChains() const									{ return mChains; }

	template<class OutputIterator>
	void				GetSequences(OutputIterator outSequences) const;

	MResidue*			GetResidue(char inChainID, uint16 inSeqNumber, char inInsertionCode);

	// statistics
	uint32				GetNrOfHBondsInParallelBridges() const				{ return mNrOfHBondsInParallelBridges; }
	uint32				GetNrOfHBondsInAntiparallelBridges() const			{ return mNrOfHBondsInAntiparallelBridges; }

	void				GetResiduesPerAlphaHelixHistogram(uint32 outHistogram[30]) const;
	void				GetParallelBridgesPerLadderHistogram(uint32 outHistogram[30]) const;
	void				GetAntiparallelBridgesPerLadderHistogram(uint32 outHistogram[30]) const;
	void				GetLaddersPerSheetHistogram(uint32 outHistogram[30]) const;
	
  private:

	void				AddResidue(const std::vector<MAtom>& inAtoms);

	void				CalculateHBondEnergies(const std::vector<MResidue*>& inResidues);
	void				CalculateAlphaHelices(const std::vector<MResidue*>& inResidues);
	void				CalculateBetaSheets(const std::vector<MResidue*>& inResidues);
	void				CalculateAccessibilities(const std::vector<MResidue*>& inResidues);

	// a thread entry point
	void				CalculateAccessibility(MResidueQueue& inQueue,
							const std::vector<MResidue*>& inResidues);

	std::string			mID, mHeader;

	std::vector<std::string>
						mCompound, mSource, mAuthor;
	std::vector<MChain*>mChains;
	uint32				mResidueCount, mChainBreaks;
	
	std::vector<std::pair<MResidue*,MResidue*> >
						mSSBonds;
	uint32				mIgnoredWaterMolecules;
	
	// statistics
	uint32				mNrOfHBondsInParallelBridges, mNrOfHBondsInAntiparallelBridges;
	uint32				mParallelBridgesPerLadderHistogram[kHistogramSize];
	uint32				mAntiparallelBridgesPerLadderHistogram[kHistogramSize];
	uint32				mLaddersPerSheetHistogram[kHistogramSize];
};

// inlines

// GetSequences can be used to quickly get all sequences in a vector<string> e.g.
template<class OutputIterator>
void MProtein::GetSequences(OutputIterator outSequences) const
{
	for (std::vector<MChain*>::const_iterator chain = mChains.begin(); chain != mChains.end(); ++chain)
	{
		std::string seq;
		(*chain)->GetSequence(seq);
		*outSequences++ = seq;
	}
}

