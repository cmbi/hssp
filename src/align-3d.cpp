// 3d dingen

//#include "mas.h"

#include <iostream>
#include <iomanip>
#include <string>
#include <limits>
#include <cmath>
#include <numeric>
#include <vector>

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
	kUnknown,
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
	kUnknown,
	
	//
	kAlanine,				// A
	kArginine,				// R
	kAsparagine,			// N
	kAsparticAcid,			// D
	kCysteine,				// C
	kGlutamicAcid,			// E
	kGlutamine,				// Q
	kGlycine,				// G
	kHistidine,				// H
	kIsoleucine,			// I
	kLeucine,				// L
	kLysine,				// K
	kMethionine,			// M
	kPhenylalanine,			// F
	kProline,				// P
	kSerine,				// S
	kThreonine,				// T
	kTryptophan,			// W
	kTyrosine,				// Y
	kValine					// V
};

struct MResidue
{
	uint32				mNumber;
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
	vector<MChain>		mChains;
};

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
			if (cAlhpaOnly and line.substr(13, 4) != " CA ")
				continue;
cerr << "calpha:" << endl << line << endl;
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
	}
	catch (const exception& e)
	{
		cerr << e.what() << endl;
		exit(1);
	}
	
	return 0;
}
