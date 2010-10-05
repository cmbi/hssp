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

#include <boost/array.hpp>

#include <boost/numeric/ublas/lu.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/lambda/lambda.hpp>

using namespace std;
using namespace tr1;

namespace fs = boost::filesystem;
namespace po = boost::program_options;
namespace io = boost::iostreams;
namespace ba = boost::algorithm;
namespace bu = boost::numeric::ublas;
namespace bl = boost::lambda;

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

/** General matrix determinant.
 * It uses lu_factorize in uBLAS. 
 */ 
template<class M> 
typename M::value_type lu_det(const M& m)
{
	// create a working copy of the input 
	M mLu(m);
	bu::permutation_matrix<uint32> pivots(m.size1());

	bu::lu_factorize(mLu, pivots);

	typename M::value_type det = 1.0;
	for (uint32 i = 0; i < pivots.size(); ++i)
	{
		if (pivots(i) != i)
			det *= -1.0;
		det *= mLu(i,i);
	} 
	return det; 
}

template<class M>
void cofactors(const M& m, M& cf)
{
	const uint32 ixs[4][3] =
	{
		{ 1, 2, 3 },
		{ 0, 2, 3 },
		{ 0, 1, 3 },
		{ 0, 1, 2 }
	};

	for (uint32 x = 0; x < 4; ++x)
	{
		const uint32* ix = ixs[x];
		
		for (uint32 y = 0; y < 4; ++y)
		{
			const uint32* iy = ixs[y];
			
			cf(x, y) =
				m(ix[0], iy[0]) * m(ix[1], iy[1]) * m(ix[2], iy[2]) +
				m(ix[0], iy[1]) * m(ix[1], iy[2]) * m(ix[2], iy[0]) +
				m(ix[0], iy[2]) * m(ix[1], iy[0]) * m(ix[2], iy[1]) -
				m(ix[0], iy[2]) * m(ix[1], iy[1]) * m(ix[2], iy[0]) -
				m(ix[0], iy[1]) * m(ix[1], iy[0]) * m(ix[2], iy[2]) -
				m(ix[0], iy[0]) * m(ix[1], iy[2]) * m(ix[2], iy[1]);
		}
	}
}

template<class M> 
typename M::value_type scale_matrix(M& m)
{
	double s = 1.0;
	
	for (uint32 x = 0; x < m.size1(); ++x)
	{
		for (uint32 y = 0; y < m.size2(); ++y)
		{
			if (s < abs(m(x, y)))
				s = abs(m(x, y));
		}
	}

	for (uint32 x = 0; x < m.size1(); ++x)
	{
		for (uint32 y = 0; y < m.size2(); ++y)
			m(x, y) /= s;
	}
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

// The next function returns the largest solution for a quartic equation
// based on Ferrari's algorithm.
// A depressed quartic is of the form:
//
//   x^4 + ax^2 + bx + c = 0
//
double largest_depressed_quartic_solution(double a, double b, double c)
{
	complex<double> P = - (a * a) / 12 - c;
	complex<double> Q = - (a * a * a) / 108 + (a * c) / 3 - (b * b) / 8;
	complex<double> R = - Q / 2.0 + sqrt((Q * Q) / 4.0 + (P * P * P) / 27.0);
	
	complex<double> U = pow(R, 1 / 3.0f);
	
	complex<double> y;
	if (U == 0.0)
		y = -5.0 * a / 6.0 + U - pow(Q, 1.0 / 3.0);
	else
		y = -5.0 * a / 6.0 + U - P / (3.0 * U);

	complex<double> W = sqrt(a + 2.0 * y);
	
	// And to get the final result:
	// result = (±W + sqrt(-(3 * alpha + 2 * y ± 2 * beta / W))) / 2;
	// We want the largest result, so:

	complex<double> t[4] =
	{
		( W + sqrt(-(3.0 * a + 2.0 * y + 2.0 * b / W))) / 2.0,
		( W + sqrt(-(3.0 * a + 2.0 * y - 2.0 * b / W))) / 2.0,
		(-W + sqrt(-(3.0 * a + 2.0 * y + 2.0 * b / W))) / 2.0,
		(-W + sqrt(-(3.0 * a + 2.0 * y - 2.0 * b / W))) / 2.0
	};
	
	cout << "0: " << t[0] << endl
		 << "1: " << t[1] << endl
		 << "2: " << t[2] << endl
		 << "3: " << t[3] << endl;

	return t[0].real();
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

	if (cAlphaA.size() != cAlphaB.size())
		throw logic_error("Protein A and B should have the same number of c-alhpa");
	
	// First calculate M, a 3x3 matrix containing the sums of products of the coordinates of A and B
	
	bu::matrix<double> M(3, 3, 0);

	for (uint32 i = 0; i < cAlphaB.size(); ++i)
	{
		point& a = cAlphaA[i];
		point& b = cAlphaB[i];
		
		M(0, 0) += a.m_x * b.m_x;
		M(0, 1) += a.m_x * b.m_y;
		M(0, 2) += a.m_x * b.m_z;
		M(1, 0) += a.m_y * b.m_x;
		M(1, 1) += a.m_y * b.m_y;
		M(1, 2) += a.m_y * b.m_z;
		M(2, 0) += a.m_z * b.m_x;
		M(2, 1) += a.m_z * b.m_y;
		M(2, 2) += a.m_z * b.m_z;
	}

//	cout << "M: " << M << endl;
//	
//	scale_matrix(M);
	
	cout << "M: " << M << endl;
	
	// Now calculate N, a 4x4 matrix
	
	bu::matrix<double> N(4, 4);
	
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

	cout << "N: " << N << endl;
	
	// det(N - λI) = 0
	// find the largest λ (λm)
	//
	// λ4 + c3λ3 + c2λ2 + c1λ + c0 = 0
	// Aλ4 + Bλ3 + Cλ2 + Dλ + E = 0
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
	
	double E = lu_det(N);
	
	cout << "C: " << C << ", D: " << D << ", E: " << E << endl;
	
	// solve quartic
	double lm = largest_depressed_quartic_solution(C, D, E);
	
	cout << "lm: " << lm << endl;
	
	bu::matrix<double> li(4, 4);
	
	li = bu::identity_matrix<double>(4) * lm;
	bu::matrix<double> t = N - li;
	
	// We're looking for em for which is said:
	// (N - lm I) em = 0
	// take em.0 = 1 so:
	
	// another approach, calculate a matrix of cofactors for t

	cout << t << endl;

	long unsigned int ixs[4][3] =
	{
		{ 1, 2, 3 },
		{ 0, 2, 3 },
		{ 0, 1, 3 },
		{ 0, 1, 2 }
	};
	
	bu::matrix<double> cf(4, 4);
	
	cofactors(t, cf);
	
	cout << "cf: " << cf << endl;
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
	
		cout << "ferrari test: " << largest_depressed_quartic_solution(6, -60, 36) << endl;
	
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
