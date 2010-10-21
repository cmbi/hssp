// A DSSP reimplementation
//
//	Copyright, M.L. Hekkelman, UMC St. Radboud, Nijmegen
//

#include "mas.h"

#include <iostream>

#include <boost/program_options.hpp>
#include <boost/filesystem/fstream.hpp>
#include <boost/foreach.hpp>
#define foreach BOOST_FOREACH
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/copy.hpp>
#include <boost/iostreams/filter/bzip2.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/date_time/gregorian/gregorian.hpp>
#include <boost/date_time/date_clock_device.hpp>

#include "structure.h"

using namespace std;
namespace fs = boost::filesystem;
namespace po = boost::program_options;
namespace io = boost::iostreams;
namespace ba = boost::algorithm;

int VERBOSE;

void ResidueToDSSPLine(const MProtein& protein, const MChain& chain,
	const MResidue& residue, ostream& os)
{
/*   
  #  RESIDUE AA STRUCTURE BP1 BP2  ACC     N-H-->O    O-->H-N    N-H-->O    O-->H-N    TCO  KAPPA ALPHA  PHI   PSI    X-CA   Y-CA   Z-CA 
 */
	boost::format kDSSPResidueLine(
	"%5.5d%5.5d%c%c %c  %c %c%c%c%c%c%c%c%4.4d%4.4d%c%4.4d %11s%11s%11s%11s  %6.3f%6.1f%6.1f%6.1f%6.1f %6.1f %6.1f %6.1f");
	
	const MAtom& ca = residue.GetCAlpha();
	
	char code = kResidueInfo[residue.GetType()].code;
	if (residue.GetType() == kCysteine and residue.GetSSBridgeNr() != 0)
		code = 'a' - 1 + (residue.GetSSBridgeNr() % 26);

	double alpha;
	char chirality;
	tr1::tie(alpha,chirality) = residue.Alpha();
	
	char ss;
	switch (residue.GetSecondaryStructure())
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
	
	string NHO[2], ONH[2];
	const HBond* acceptors = residue.Acceptor();
	const HBond* donors = residue.Donor();
	for (uint32 i = 0; i < 2; ++i)
	{
		NHO[i] = ONH[i] = "0, 0.0";
		
		if (acceptors[i].residue != nil)
		{
			int32 d = acceptors[i].residue->GetNumber() - residue.GetNumber();
			NHO[i] = (boost::format("%d,%3.1f") % d % acceptors[i].energy).str();
		}
	
		if (donors[i].residue != nil)
		{
			int32 d = donors[i].residue->GetNumber() - residue.GetNumber();
			ONH[i] = (boost::format("%d,%3.1f") % d % donors[i].energy).str();
		}
	}
	
	uint32 bp[2] = {};
	char bridgelabel[2] = { ' ', ' ' };
	for (uint32 i = 0; i < 2; ++i)
	{
		MBridgeParner p = residue.GetBetaPartner(i);
		if (p.residue != nil)
		{
			bp[i] = p.residue->GetNumber();
			bridgelabel[i] = 'A' + p.ladder % 26;
			if (p.parallel)
				bridgelabel[i] = tolower(bridgelabel[i]);
		}
	}
	
	char sheet = ' ';
	if (residue.GetSheet() != 0)
		sheet = 'A' + (residue.GetSheet() - 1) % 26;
	
	char helix[3];
	for (uint32 stride = 3; stride <= 5; ++stride)
	{
		switch (residue.GetHelixFlag(stride))
		{
			case helixNone:			helix[stride - 3] = ' '; break;
			case helixStart:		helix[stride - 3] = '>'; break;
			case helixEnd:			helix[stride - 3] = '<'; break;
			case helixStartAndEnd:	helix[stride - 3] = 'X'; break;
			case helixMiddle:		helix[stride - 3] = '0' + stride; break;
		}
	}
	
	char bend = ' ';
	if (residue.IsBend())
		bend = 'S';
	
	os << (kDSSPResidueLine % residue.GetNumber() % ca.mResSeq % ca.mICode % ca.mChainID % code %
		ss % helix[0] % helix[1] % helix[2] % bend % chirality % bridgelabel[0] % bridgelabel[1] %
		bp[0] % bp[1] % sheet % residue.Accessibility() %
		NHO[0] % ONH[0] % NHO[1] % ONH[1] %
		residue.TCO() % residue.Kappa() % alpha % residue.Phi() % residue.Psi() %
		ca.mLoc.mX % ca.mLoc.mY % ca.mLoc.mZ) << endl;
}

void ReportDSSP(MProtein& protein, ostream& os)
{
	const string kFirstLine("==== Secondary Structure Definition by the program DSSP, CMBI version by M.L. Hekkelman/2010-10-21 ==== ");
	boost::format kHeaderLine("%1% %|127t|%2%");
	
	using namespace boost::gregorian;
	
	uint32 nrOfResidues, nrOfChains, nrOfSSBridges, nrOfIntraChainSSBridges;
	protein.GetStatistics(nrOfResidues, nrOfChains, nrOfSSBridges, nrOfIntraChainSSBridges);
	
	date today = day_clock::local_day();

	os << kHeaderLine % (kFirstLine + "DATE=" + to_iso_extended_string(today)) % '.' << endl;
	os << kHeaderLine % "REFERENCE W. KABSCH AND C.SANDER, BIOPOLYMERS 22 (1983) 2577-2637" % '.' << endl;
	os << kHeaderLine % protein.GetHeader() % '.' << endl;
	os << kHeaderLine % protein.GetCompound() % '.' << endl;
	os << kHeaderLine % protein.GetSource() % '.' << endl;
	os << kHeaderLine % protein.GetAuthor() % '.' << endl;

	os << boost::format("%5.5d%3.3d%3.3d%3.3d%3.3d TOTAL NUMBER OF RESIDUES, NUMBER OF CHAINS, NUMBER OF SS-BRIDGES(TOTAL,INTRACHAIN,INTERCHAIN) %|127t|%c") %
   		nrOfResidues % nrOfChains % nrOfSSBridges % nrOfIntraChainSSBridges % (nrOfSSBridges - nrOfIntraChainSSBridges) % '.' << endl;
//	   << (boost::format("%8.1f   ACCESSIBLE SURFACE OF PROTEIN (ANGSTROM**2)") %
//	   		protein.GetAccessibleSurface()) << endl;

	os << "  #  RESIDUE AA STRUCTURE BP1 BP2  ACC     N-H-->O    O-->H-N    N-H-->O    O-->H-N    TCO  KAPPA ALPHA  PHI   PSI    X-CA   Y-CA   Z-CA " << endl;
	boost::format kDSSPResidueLine(
		"%5.5d        !*             0   0    0      0, 0.0     0, 0.0     0, 0.0     0, 0.0   0.000 360.0 360.0 360.0 360.0    0.0    0.0    0.0");

	foreach (const MChain* chain, protein.GetChains())
	{
		foreach (const MResidue* residue, chain->GetResidues())
			ResidueToDSSPLine(protein, *chain, *residue, os);

		if (chain != protein.GetChains().back())
			os << (kDSSPResidueLine % (chain->GetResidues().back()->GetNumber() + 1)) << endl;
	}
}

int main(int argc, char* argv[])
{
	try
	{
		po::options_description desc("DSSP options");
		desc.add_options()
			("help,h",							 "Display help message")
			("input,i",		po::value<string>(), "Input file")
			("output,o",	po::value<string>(), "Output file, use 'stdout' to output to screen")
			("verbose,v",						 "Verbose output")
			;
	
		po::positional_options_description p;
		p.add("input", 1);
		p.add("output", 2);
	
		po::variables_map vm;
		po::store(po::command_line_parser(argc, argv).options(desc).positional(p).run(), vm);
		po::notify(vm);

		if (vm.count("help") or not vm.count("input"))
		{
			cerr << desc << endl;
			exit(1);
		}

		VERBOSE = vm.count("verbose");
		if (vm.count("debug"))
			VERBOSE = vm["debug"].as<int>();
		
		string input = vm["input"].as<string>();

		ifstream infile(input.c_str(), ios_base::in | ios_base::binary);
		if (not infile.is_open())
			throw runtime_error("No such file");
		
		io::filtering_stream<io::input> in;
		
		if (ba::ends_with(input, ".bz2"))
			in.push(io::bzip2_decompressor());
		else if (ba::ends_with(input, ".gz"))
			in.push(io::gzip_decompressor());
		
		in.push(infile);
	
		// OK, we've got the file, now create a protein
		MProtein a(in);
		
		// then calculate the secondary structure
		a.CalculateSecondaryStructure();
	
		// and finally report the secondary structure in the DSSP format
		// either to cout or an (optionally compressed) file.
		if (vm.count("output"))
		{
			string output = vm["output"].as<string>();
			
			ofstream outfile(output.c_str(), ios_base::out|ios_base::trunc|ios_base::binary);
			if (not outfile.is_open())
				throw runtime_error("could not create output file");
			
			io::filtering_stream<io::output> out;
			if (ba::ends_with(output, ".bz2"))
				out.push(io::bzip2_compressor());
			else if (ba::ends_with(output, ".gz"))
				out.push(io::gzip_compressor());
			
			out.push(outfile);
			
			ReportDSSP(a, out);
		}
		else
			ReportDSSP(a, cout);
	}
	catch (const exception& e)
	{
		cerr << e.what() << endl;
		exit(1);
	}
	
	return 0;
}

