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

#include <boost/numeric/ublas/lu.hpp>
#include <boost/numeric/ublas/io.hpp>

using namespace std;
using namespace tr1;
namespace fs = boost::filesystem;
namespace po = boost::program_options;
namespace io = boost::iostreams;
namespace ba = boost::algorithm;
namespace bu = boost::numeric::ublas;

// --------------------------------------------------------------------

struct point
{

				point(const point& rhs)
					: m_x(rhs.m_x), m_y(rhs.m_y), m_z(rhs.m_z) {}

				point() : m_x(0), m_y(0), m_z(0) {}

	point&		operator=(const point& rhs)
				{
					m_x = rhs.m_x;
					m_y = rhs.m_y;
					m_z = rhs.m_z;
					
					return *this;
				}

	float		m_x, m_y, m_z;
};

//      /** General matrix inversion routine. 
//       * It uses lu_factorize and lu_substitute in uBLAS to invert a matrix 
//       */ 
//      template<class M1, class M2> 
//      void lu_inv(M1 const& m, M2& inv) { 
//        JFR_PRECOND(m.size1() == m.size2(), 
//                    "ublasExtra::lu_inv(): input matrix must be squared"); 
//        JFR_PRECOND(inv.size1() == m.size1() && inv.size1() == m.size2(), 
//                    "ublasExtra::lu_inv(): invalid size for inverse matrix"); 
//
//        using namespace boost::numeric::ublas; 
//        // create a working copy of the input 
//        mat mLu(m); 
//
//        // perform LU-factorization 
//        lu_factorize(mLu); 
//
//        // create identity matrix of "inverse" 
//        inv.assign(identity_mat(m.size1())); 
//
//        // backsubstitute to get the inverse 
//        lu_substitute<mat const, M2 >(mLu, inv); 
//      } 

/** General matrix determinant. 
 * It uses lu_factorize in uBLAS. 
 */ 
template<class M> 
float lu_det(const M& m)
{
	// create a working copy of the input 
	bu::matrix<float> mLu(m);
	bu::permutation_matrix<uint32> pivots(m.size1());

	bu::lu_factorize(mLu, pivots);

	float det = 1.0;
	for (uint32 i = 0; i < pivots.size(); ++i)
	{
		if (pivots(i) != i)
			det *= -1.0;
		det *= mLu(i,i);
	} 
	return det; 
}

//template<class M>
//float cofactor(const M& m, uint32 x, uint32 y)
//{
//	
//}

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
};

struct MAtom
{
	MAtomType			mType;
	point				mLoc;
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
	
	void				GetCAlphaLocations(char inChain, vector<point>& outPoints) const;
};

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
			if (a.mType == kCarbon)
				outPoints.push_back(a.mLoc);
		}
	}
}

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
				os << "Atom: " << atom.mLoc.m_x << ", " << atom.mLoc.m_y << ", " << atom.mLoc.m_z << endl;
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
			continue;
		}
		
		if (ba::starts_with(line, "ATOM  "))
		{
			if (cAlhpaOnly and line.substr(12, 4) != " CA ")
				continue;
			
			MAtom atom;
			atom.mType = kCarbon;
			atom.mLoc.m_x = ParseFloat(line.substr(30, 8));
			atom.mLoc.m_y = ParseFloat(line.substr(38, 8));
			atom.mLoc.m_z = ParseFloat(line.substr(46, 8));
			
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

float largest_quartic_solution(float A, float B, float C, float D, float E)
{
	// See http://en.wikipedia.org/wiki/Quartic_function#Solving_a_quartic_equation
	float alpha = -3 * (B * B) / 8 * A * A + C / A;
//	float alpha = C;		// since B == 0 and A == 1.0
	
	float beta = B * B * B / 8 * A * A * A - B * C / 2 * A * A + D / A;
//	float beta = D;
	
	float gamma = -3 * B * B * B * B / 256 * A * A * A * A + C * B * B / 16 * A * A * A -
					B * D / 4 * A * A + E / A;
//	float gamma = E;
	cout << "alpha: " << alpha << ", beta: " << beta << ", gamma: " << gamma << endl;

	float P = - (alpha * alpha) / 12 - gamma;
	float Q = - (alpha * alpha * alpha) / 108 + (alpha * gamma) / 3 - (beta * beta) / 8;
	float R = - Q / 2 + sqrt((Q * Q) / 4 + (P * P * P) / 27);
	
	float U = pow(R, 1 / 3.0f);
	
	cout << "P: " << P << ", Q: " << Q << ", R: " << R << ", U: " << U << endl;
	
	float y;
	if (U == 0)
		y = -5 * alpha / 6 + U - pow(Q, 1 / 3.0f);
	else
		y = -5 * alpha / 6 + U - P / (3 * U);

	cout << "y: " << y << endl;
	
	float W = sqrt(alpha + 2 * y);
	cout << "W: " << W << endl;
	
	float result;

	float x = (W + sqrt(-(3 * alpha + 2 * y + 2 * beta / W))) / 2;
	result = x;

	x = (W - sqrt(-(3 * alpha + 2 * y + 2 * beta / W))) / 2;
	if (result < x or isnan(result))
		result = x;

	x = (-W + sqrt(-(3 * alpha + 2 * y - 2 * beta / W))) / 2;
	if (result < x or isnan(result))
		result = x;

	x = (-W - sqrt(-(3 * alpha + 2 * y - 2 * beta / W))) / 2;
	if (result < x or isnan(result))
		result = x;
}

void align_proteins(MProtein& a, MProtein& b)
{
	vector<point> cAlphaA;
	a.GetCAlphaLocations(0, cAlphaA);
	
	point translationA = center_points(cAlphaA);
	cout << "translation to centroid for A: " << translationA.m_x << "," << translationA.m_y << "," << translationA.m_z << endl;
	
	vector<point> cAlphaB;
	b.GetCAlphaLocations(0, cAlphaB);
	
	point translationB = center_points(cAlphaB);
	cout << "translation to centroid for B: " << translationB.m_x << "," << translationB.m_y << "," << translationB.m_z << endl;
	
	// First calculate M, a 3x3 matrix containing the sums of products of the coordinates of A and B
	
	bu::matrix<float> M(3, 3);
	
	foreach (point& a, cAlphaA)
	{
		foreach (point& b, cAlphaB)
		{
			M(0, 0) = a.m_x * b.m_x;
			M(0, 1) = a.m_x * b.m_y;
			M(0, 2) = a.m_x * b.m_z;
			M(1, 0) = a.m_y * b.m_x;
			M(1, 1) = a.m_y * b.m_y;
			M(1, 2) = a.m_y * b.m_z;
			M(2, 0) = a.m_z * b.m_x;
			M(2, 1) = a.m_z * b.m_y;
			M(2, 2) = a.m_z * b.m_z;
		}
	}
	
	cout << "M: " << endl << M << endl << endl;

	// Now calculate N, a 4x4 matrix
	
	bu::matrix<float> N(4, 4);
	
	N(0, 0)				= M(0, 0) + M(1, 1) + M(2, 2);
	N(0, 1)	= N(1, 0)	= M(1, 2) - M(2, 1);
	N(0, 2) = N(2, 0)	= M(2, 0) - M(0, 2);
	N(0, 3) = N(3, 0)	= M(0, 1) - M(1, 0);
	
	N(1, 1)				= M(0, 0) - M(1, 1) - M(2, 2);
	N(1, 2) = N(2, 1)	= M(0, 1) + M(1, 0);
	N(1, 3) = N(3, 1)	= M(0, 2) + M(2, 0);
	
	N(2, 2)				= -M(0, 0) + M(1, 1) - M(2, 2);
	N(2, 3) = N(3, 2)	= M(1, 2) + M(2, 1);
	
	N(3, 3)				= -M(0, 0) - M(1, 1) + M(2, 2);

	cout << "N: " << endl << N << endl << endl;
	
	// det(N - λI) = 0
	// find the largest λ (λm)
	//
	// λ4 + c3λ3 + c2λ2 + c1λ + c0 = 0
	// Aλ4 + Bλ3 + Cλ2 + Dλ + E = 0
	// B = 0
	// and so this is a so-called depressed quartic
	// solve it using Ferrari's algorithm
	
	float A = 1;
	
	float B = 0;
	
	float C = -2 * (
		M(0, 0) * M(0, 0) + M(0, 1) * M(0, 1) + M(0, 2) * M(0, 2) +
		M(1, 0) * M(1, 0) + M(1, 1) * M(1, 1) + M(1, 2) * M(1, 2) +
		M(2, 0) * M(2, 0) + M(2, 1) * M(2, 1) + M(2, 2) * M(2, 2));
	
	float D = 8 * (M(0, 0) * M(1, 2) * M(2, 1) +
					M(1, 1) * M(2, 0) * M(0, 2) +
					M(2, 2) * M(0, 1) * M(1, 0)) -
			   8 * (M(0, 0) * M(1, 1) * M(2, 2) +
					M(1, 2) * M(2, 0) * M(0, 1) +
					M(2, 1) * M(1, 0) * M(0, 2));
	
	float E = lu_det(N);
	
	// solve quartic
	
	float lm = largest_quartic_solution(A, B, C, D, E);
	
	bu::matrix<float> li(4, 4);
	
	li = bu::identity_matrix<float>(4) * lm;
	bu::matrix<float> t = N - li;
	
	// We're looking for em for which is said:
	// (N - lm I) em = 0
	// take em.0 = 1 so:
	
	
}

float solve_cubic(float A, float B, float C)
{
	
}

int main(int argc, char* argv[])
{
	try
	{
		po::options_description desc("align-3d options");
		desc.add_options()
			("help,h",							 			"Display help message")
			("input,i",		po::value<vector<string> >(),	"Input files")
			;
	
		po::positional_options_description p;
		p.add("input", -1);
	
		po::variables_map vm;
		po::store(po::command_line_parser(argc, argv).options(desc).positional(p).run(), vm);
		po::notify(vm);
		
		if (vm.count("help") or vm.count("input") == 0 or vm["input"].as<vector<string> >().size() != 2)
		{
			cerr << desc << endl;
			exit(1);
		}
		
		MProtein a, b;
		
		fs::ifstream file_a(vm["input"].as<vector<string> >()[0]);
		ParsePDB(file_a, a, true);

		fs::ifstream file_b(vm["input"].as<vector<string> >()[1]);
		ParsePDB(file_b, b, true);
		
		align_proteins(a, b);
	}
	catch (const exception& e)
	{
		cerr << e.what() << endl;
		exit(1);
	}
	
	return 0;
}
