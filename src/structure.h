// structure stuff for proteins

#pragma once

#include "primitives-3d.h"

struct MAtom;
class MResidue;
class MChain;
class MProtein;

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

MAtomType MapElement(const std::string& inElement);

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
	void		Rotate(const MQuaternion& inRotation)	{ mLoc.rotate(inRotation); }
	void		WritePDB(std::ostream& os);

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

extern const MResidueInfo kResidueInfo[];

char MapThreeLetterCode(const std::string& inThreeLetterCode);

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

enum MSecondaryStructure
{
	loop, alphahelix, betabridge, strand, helix_3, helix_5, turn, bend
};

class MResidue
{
  public:
						MResidue(MChain& chain, uint32 inSeqNumber, uint32 inNumber,
							MResidueType inType,
							MResidue* inPrevious, const std::vector<MAtom>& inAtoms);

	void				SetChainID(char inID);
	char				GetChainID() const;

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
	
	void				SetSecondaryStructure(MSecondaryStructure inSS)
														{ mSecondaryStructure = inSS; }
	MSecondaryStructure	GetSecondaryStructure() const	{ return mSecondaryStructure; }
	
	const MResidue*		Next() const					{ return mPrev; }
	const MResidue*		Prev() const					{ return mNext; }
	
	void				SetBetaPartner(uint32 n, MResidue* inResidue, uint32 inLadder)
													{ mBetaPartner[n] = inResidue; }
	const MResidue*		GetBetaPartner(uint32 n) const
													{ return mBetaPartner[n]; }
	void				SetSheet(uint32 inSheet)	{ mSheet = inSheet; }

	void				SetSSBridgeNr(uint8 inBridgeNr);
	uint8				GetSSBridgeNr() const;

	void				AddAtom(MAtom& inAtom);
	
	HBond*				Donor()						{ return mHBondDonor; }
	HBond*				Acceptor()					{ return mHBondAcceptor; }

	bool				ValidDistance(const MResidue& inNext) const;

	static bool			TestBond(const MResidue* a, const MResidue* b)
						{
							return a->TestBond(b);
						}

	// bridge functions
	MBridgeType			TestBridge(MResidue* inResidue) const;

	uint16				GetSeqNumber() const		{ return mSeqNumber; }
	uint16				GetNumber() const			{ return mNumber; }

	void				Translate(const MPoint& inTranslation);
	void				Rotate(const MQuaternion& inRotation);

	void				WritePDB(std::ostream& os);
	void				WriteDSSP(std::ostream& os);

	static double		CalculateHBondEnergy(MResidue& inDonor, MResidue& inAcceptor);

	std::vector<MAtom>&	GetAtoms()					{ return mAtoms; }
	const std::vector<MAtom>&
						GetAtoms() const			{ return mAtoms; }

  protected:

	bool				TestBond(const MResidue* other) const;

	MChain&				mChain;
	MResidue*			mPrev;
	MResidue*			mNext;
	int32				mSeqNumber, mNumber;
	MResidueType		mType;
	uint8				mSSBridgeNr;
	MSecondaryStructure	mSecondaryStructure;
	MAtom				mC, mN, mCA, mO, mH;
	HBond				mHBondDonor[2], mHBondAcceptor[2];
	std::vector<MAtom>	mAtoms;
	MResidue*			mBetaPartner[2];
	uint32				mSheet;
};

class MChain
{
  public:

						MChain(char inChainID = 0) : mChainID(inChainID) {}

	char				GetChainID() const					{ return mChainID; }
	void				SetChainID(char inID);

	MResidue&			GetResidueBySeqNumber(uint16 inSeqNumber);

	void				Translate(const MPoint& inTranslation);
	void				Rotate(const MQuaternion& inRotation);

	void				WritePDB(std::ostream& os);
	void				WriteDSSP(std::ostream& os);
	
	std::vector<MResidue*>&
						GetResidues()						{ return mResidues; }
	const std::vector<MResidue*>&
						GetResidues() const					{ return mResidues; }

	void				AddResidue(uint32 inSeqNumber, uint32 inNumber,
							MResidueType inType, const std::vector<MAtom>& inAtoms);

  private:
	char				mChainID;
	std::vector<MResidue*>
						mResidues;
};

class MProtein
{
  public:
						MProtein() {}

	std::string			GetID() const					{ return mID; }
						
						MProtein(std::istream& is, bool inCAlphaOnly = false);

	void				CalculateSecondaryStructure();
	void				CalculateSSBridges();
	
	void				GetCAlphaLocations(char inChain, std::vector<MPoint>& outPoints) const;
	MPoint				GetCAlphaPosition(char inChain, int16 inPDBResSeq) const;
	
	void				GetSequence(char inChain, entry& outEntry) const;
	void				GetSequence(char inChain, sequence& outSequence) const;

	void				Center();
	void				Translate(const MPoint& inTranslation);
	void				Rotate(const MQuaternion& inRotation);

	void				WritePDB(std::ostream& os);

	void				WriteDSSP(std::ostream& os);
	
	void				GetPoints(std::vector<MPoint>& outPoints) const;

	MResidue&			GetResidue(MResidueID inID)							{ return GetResidue(inID.chain, inID.seqNumber); }
	MResidue&			GetResidue(char inChainID, uint16 inSeqNumber);

	MChain&				GetChain(char inChainID);

	std::map<char,MChain>&
						GetChains() 						{ return mChains; }
	
  private:

	std::string			mID;
	std::map<char,MChain>
						mChains;
	
	std::vector<std::pair<MResidueID,MResidueID>>
						mSSBonds;
};
