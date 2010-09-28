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

using namespace std;
using namespace tr1;
namespace fs = boost::filesystem;
namespace po = boost::program_options;
namespace io = boost::iostreams;
namespace ba = boost::algorithm;

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
};

struct MAtom
{
	MAtomType			mType;
	float				mX, mY, mZ;
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
};

struct MChain
{
	char				mChainID;
	vector<MResidue>	mResidues;
};

struct MProtein
{
	string				mID;
	map<char,MChain>	mChains;
};

ostream& operator<<(ostream& os, const MProtein& prot)
{
	os << "ID: " << prot.mID << endl;
	
	for (map<char,MChain>::const_iterator c = prot.mChains.begin(); c != prot.mChains.end(); ++c)
	{
		const MChain& chain = c->second;
		os << "Chain: " << chain.mChainID << endl << endl;
		
		foreach (const MResidue& res, chain.mResidues)
		{
			os << "Residue: " << res.mNumber << ' ' << kResidueInfo[res.mType].code << ' ' << kResidueInfo[res.mType].name << endl;
			
			foreach (const MAtom& atom, res.mAtoms)
			{
				os << "Atom: " << atom.mX << ", " << atom.mY << ", " << atom.mZ << endl;
			}
		}
	}

	return os;	
}

inline float ParseFloat(const string& s)
{
	return boost::lexical_cast<float>(ba::trim_copy(s));
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
		}
		
		if (ba::starts_with(line, "ATOM  "))
		{
			if (cAlhpaOnly and line.substr(12, 4) != " CA ")
				continue;
			
			MAtom atom;
			atom.mType = kCarbon;
			atom.mX = ParseFloat(line.substr(30, 8));
			atom.mY = ParseFloat(line.substr(38, 8));
			atom.mZ = ParseFloat(line.substr(46, 8));
			
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

int main(int argc, char* argv[])
{
	try
	{
		po::options_description desc("align-3d options");
		desc.add_options()
			("help,h",							 "Display help message")
			("input,i",		po::value<string>(), "Input file")
			;
	
		po::positional_options_description p;
		p.add("input", -1);
	
		po::variables_map vm;
		po::store(po::command_line_parser(argc, argv).options(desc).positional(p).run(), vm);
		po::notify(vm);
		
		if (vm.count("help") or (vm.count("input") == 0))
		{
			cerr << desc << endl;
			exit(1);
		}
		
		fs::ifstream file(vm["input"].as<string>());
		MProtein prot;
		
		ParsePDB(file, prot, true);
		
		cout << prot << endl;
	}
	catch (const exception& e)
	{
		cerr << e.what() << endl;
		exit(1);
	}
	
	return 0;
}
