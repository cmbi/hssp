// Copyright Maarten L. Hekkelman, Radboud University 2008-2011.
// Copyright Coos Baakman, Jon Black, Wouter G. Touw & Gert Vriend, Radboud university medical center 2015.
//   Distributed under the Boost Software License, Version 1.0.
//       (See accompanying file LICENSE_1_0.txt or copy at
//             http://www.boost.org/LICENSE_1_0.txt)

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "mas.h"

#include "dssp.h"
#include "structure.h"
#include "utils.h"

#if defined(_MSC_VER)
#include <conio.h>
#include <ctype.h>
#endif
#include <iostream>
#include <ctime>
#include <algorithm>


std::string LimitWidth(int64 n, uint64 w)
{
    std::string s = std::to_string(n);

    if (s.length() > w)
    {
        s = s.substr(s.length() - w);

        if (n < 0)
            s = s.replace(0, 1, "-");
    }
    else if (s.length() < w)
    {
        s = s.insert(0, w - s.length(), ' ');
    }

    return s;
}


std::string ResidueToDSSPLine(const MResidue& residue)
{
/*
  This is the header line for the residue lines in a DSSP file:

  #  RESIDUE AA STRUCTURE BP1 BP2  ACC     N-H-->O    O-->H-N    N-H-->O    O-->H-N    TCO  KAPPA ALPHA  PHI   PSI    X-CA   Y-CA   Z-CA           CHAIN AUTHCHAIN
 */
  std::string kDSSPResidueLine(
  "%5s%5s%1.1s%1.1s %c  %c %c%c%c%c%c%c%c%4s%4s%c%4.4d %11s%11s%11s%11s  %6.3f%6.1f%6.1f%6.1f%6.1f %6.1f %6.1f %6.1f             %4.4s      %4.4s");

  const MAtom& ca = residue.GetCAlpha();

  char code = kResidueInfo[residue.GetType()].code;
  if (residue.GetType() == kCysteine and residue.GetSSBridgeNr() != 0)
    code = 'a' + ((residue.GetSSBridgeNr() - 1) % 26);

  char ss = ' ';
  switch (residue.GetSecondaryStructure())
  {
    case alphahelix:  ss = 'H'; break;
    case betabridge:  ss = 'B'; break;
    case strand:    ss = 'E'; break;
    case helix_3:    ss = 'G'; break;
    case helix_5:    ss = 'I'; break;
    case turn:      ss = 'T'; break;
    case bend:      ss = 'S'; break;
    case loop:      ss = ' '; break;
  }

  char helix[3];
  for (uint32 stride = 3; stride <= 5; ++stride)
  {
    switch (residue.GetHelixFlag(stride))
    {
      case helixNone:      helix[stride - 3] = ' '; break;
      case helixStart:    helix[stride - 3] = '>'; break;
      case helixEnd:      helix[stride - 3] = '<'; break;
      case helixStartAndEnd:  helix[stride - 3] = 'X'; break;
      case helixMiddle:    helix[stride - 3] = '0' + stride; break;
    }
  }

  char bend = ' ';
  if (residue.IsBend())
    bend = 'S';

  double alpha;
  char chirality;
  std::tie(alpha,chirality) = residue.Alpha();

  std::string bp[2] = {};
  char bridgelabel[2] = { ' ', ' ' };
  for (uint32 i = 0; i < 2; ++i)
  {
    MBridgeParner p = residue.GetBetaPartner(i);
    if (p.residue != nullptr)
    {
      bp[i] = LimitWidth(p.residue->GetNumber(), 4);  // won't fit otherwise...
      bridgelabel[i] = 'A' + p.ladder % 26;
      if (p.parallel)
        bridgelabel[i] = tolower(bridgelabel[i]);
    }
  }

  char sheet = ' ';
  if (residue.GetSheet() != 0)
    sheet = 'A' + (residue.GetSheet() - 1) % 26;

  std::string NHO[2], ONH[2];
  const HBond* acceptors = residue.Acceptor();
  const HBond* donors = residue.Donor();
  for (uint32 i = 0; i < 2; ++i)
  {
    NHO[i] = ONH[i] = "0, 0.0";

    if (acceptors[i].residue != nullptr)
    {
      std::string d = LimitWidth(acceptors[i].residue->GetNumber() - residue.GetNumber(), 5);  // won't fit otherwise...
      NHO[i] = Format("%s,%3.1f", d, acceptors[i].energy);
    }

    if (donors[i].residue != nullptr)
    {
      std::string d = LimitWidth(donors[i].residue->GetNumber() - residue.GetNumber(), 5);  // won't fit otherwise...
      ONH[i] = Format("%s,%3.1f", d, donors[i].energy);
    }
  }

  std::string chainChar = ca.mChainID,
                          long_ChainID1 = ca.mChainID,
                          long_ChainID2 = ca.mAuthChainID;
  if (ca.mChainID.length () > 1)
  {
    // For mmCIF compatibility

    chainChar = ">";
  }

  // won't fit otherwise...
  std::string residueNumber = LimitWidth(residue.GetNumber(), 5),
              seqNumber = LimitWidth(ca.mResSeq, 5);

  return Format(kDSSPResidueLine, residueNumber, seqNumber, ca.mICode, chainChar, code,
                ss, helix[0], helix[1], helix[2], bend, chirality, bridgelabel[0], bridgelabel[1],
                bp[0], bp[1], sheet, floor(residue.Accessibility() + 0.5),
                NHO[0], ONH[0], NHO[1], ONH[1],
                residue.TCO(), residue.Kappa(), alpha, residue.Phi(), residue.Psi(),
                ca.mLoc.mX, ca.mLoc.mY, ca.mLoc.mZ, long_ChainID1, long_ChainID2);
}

void WriteDSSP(MProtein& protein, std::ostream& os)
{
  const std::string kFirstLine("==== Secondary Structure Definition by the program DSSP, CMBI version " PACKAGE_VERSION "                          ==== ");
  std::string kHeaderLine("%1% %|127t|%2%");

  uint32 nrOfResidues, nrOfChains, nrOfSSBridges, nrOfIntraChainSSBridges, nrOfHBonds;
  uint32 nrOfHBondsPerDistance[11] = {};

  protein.GetStatistics(nrOfResidues, nrOfChains, nrOfSSBridges, nrOfIntraChainSSBridges, nrOfHBonds, nrOfHBondsPerDistance);

  time_t today = time(0);
  char todayString[sizeof "2011-10-08T07:07:09Z"];
  strftime(todayString, sizeof(todayString), "%FT%TZ", gmtime(&today));

  os << Format(kHeaderLine, kFirstLine + "DATE=" + todayString, '.') << std::endl;
  os << Format(kHeaderLine, "REFERENCE W. KABSCH AND C.SANDER, BIOPOLYMERS 22 (1983) 2577-2637", '.') << std::endl;
  os << Format(kHeaderLine, protein.GetHeader(), '.') << std::endl;
  if (not protein.GetCompound().empty())
    os << Format(kHeaderLine, protein.GetCompound(), '.') << std::endl;
  if (not protein.GetSource().empty())
    os << Format(kHeaderLine, protein.GetSource(), '.') << std::endl;
  if (not protein.GetAuthor().empty())
    os << Format(kHeaderLine, protein.GetAuthor(), '.') << std::endl;

  double accessibleSurface = 0;  // calculate accessibility as
  for (const MChain* chain : protein.GetChains())
  {
    for (const MResidue* residue : chain->GetResidues())
      accessibleSurface += residue->Accessibility();
  }

  os << Format("%5.5d%3.3d%3.3d%3.3d%3.3d TOTAL NUMBER OF RESIDUES, NUMBER OF CHAINS, NUMBER OF SS-BRIDGES(TOTAL,INTRACHAIN,INTERCHAIN) %|127t|%c",
       nrOfResidues, nrOfChains, nrOfSSBridges, nrOfIntraChainSSBridges, (nrOfSSBridges - nrOfIntraChainSSBridges), '.') << std::endl;
     os << Format(kHeaderLine, Format("%8.1f   ACCESSIBLE SURFACE OF PROTEIN (ANGSTROM**2)", accessibleSurface), '.') << std::endl;

  // hydrogenbond summary

  os << Format(kHeaderLine, Format("%5.5d%5.1f   TOTAL NUMBER OF HYDROGEN BONDS OF TYPE O(I)-->H-N(J)  , SAME NUMBER PER 100 RESIDUES",
                                   nrOfHBonds, (nrOfHBonds * 100.0 / nrOfResidues)), '.') << std::endl;

  uint32 nrOfHBondsInParallelBridges = protein.GetNrOfHBondsInParallelBridges();
  os << Format(kHeaderLine,
    Format("%5.5d%5.1f   TOTAL NUMBER OF HYDROGEN BONDS IN     PARALLEL BRIDGES, SAME NUMBER PER 100 RESIDUES",
           nrOfHBondsInParallelBridges, nrOfHBondsInParallelBridges * 100.0 / nrOfResidues), '.') << std::endl;

  uint32 nrOfHBondsInAntiparallelBridges = protein.GetNrOfHBondsInAntiparallelBridges();
  os << Format(kHeaderLine,
    Format("%5.5d%5.1f   TOTAL NUMBER OF HYDROGEN BONDS IN ANTIPARALLEL BRIDGES, SAME NUMBER PER 100 RESIDUES",
           nrOfHBondsInAntiparallelBridges, nrOfHBondsInAntiparallelBridges * 100.0 / nrOfResidues), '.') << std::endl;

  std::string kHBondsLine("%5.5d%5.1f   TOTAL NUMBER OF HYDROGEN BONDS OF TYPE O(I)-->H-N(I%c%1.1d), SAME NUMBER PER 100 RESIDUES");
  for (int32 k = 0; k < 11; ++k)
  {
    os << Format(kHeaderLine, Format(kHBondsLine, nrOfHBondsPerDistance[k], (nrOfHBondsPerDistance[k] * 100.0 / nrOfResidues), (k - 5 < 0 ? '-' : '+'), abs(k - 5)), '.') << std::endl;
  }

  // histograms...

  uint32 histogram[kHistogramSize];
  os << "  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30     *** HISTOGRAMS OF ***           ." << std::endl;

  protein.GetResiduesPerAlphaHelixHistogram(histogram);
  for (uint32 i = 0; i < kHistogramSize; ++i)
    os << Format("%3.3d", histogram[i]);
  os << "    RESIDUES PER ALPHA HELIX         ." << std::endl;

  protein.GetParallelBridgesPerLadderHistogram(histogram);
  for (uint32 i = 0; i < kHistogramSize; ++i)
    os << Format("%3.3d", histogram[i]);
  os << "    PARALLEL BRIDGES PER LADDER      ." << std::endl;

  protein.GetAntiparallelBridgesPerLadderHistogram(histogram);
  for (uint32 i = 0; i < kHistogramSize; ++i)
    os << Format("%3.3d", histogram[i]);
  os << "    ANTIPARALLEL BRIDGES PER LADDER  ." << std::endl;

  protein.GetLaddersPerSheetHistogram(histogram);
  for (uint32 i = 0; i < kHistogramSize; ++i)
    os << Format("%3.3d", histogram[i]);
  os << "    LADDERS PER SHEET                ." << std::endl;

  // per residue information

  os << "  #  RESIDUE AA STRUCTURE BP1 BP2  ACC     N-H-->O    O-->H-N    N-H-->O    O-->H-N    TCO  KAPPA ALPHA  PHI   PSI    X-CA   Y-CA   Z-CA            CHAIN AUTHCHAIN" << std::endl;
  std::string kDSSPResidueLine("%5s        !%c             0   0    0      0, 0.0     0, 0.0     0, 0.0     0, 0.0   0.000 360.0 360.0 360.0 360.0    0.0    0.0    0.0");

  std::vector<const MResidue*> residues;

  for (const MChain* chain : protein.GetChains())
  {
    for (const MResidue* residue : chain->GetResidues())
    {
      residues.push_back(residue);
    }
  }

  // keep residues sorted by residue number as assigned during reading the PDB file
  std::sort(residues.begin(), residues.end(),
       [](const MResidue *r1, const MResidue *r2)
       {
            return r1->GetNumber() > r2->GetNumber();
       }
  );

  const MResidue* last = nullptr;
  for (const MResidue* residue : residues)
  {
    // insert a break line whenever we detect missing residues
    // can be the transition to a different chain, or missing residues in the current chain
    if (last != nullptr and last->GetNumber() + 1 != residue->GetNumber())
    {
      char breaktype = ' ';
      if (last->GetChainID() != residue->GetChainID())
        breaktype = '*';

      std::string residueNumber = LimitWidth(last->GetNumber() + 1, 5);
      os << Format(kDSSPResidueLine, residueNumber, breaktype) << std::endl;
    }
    os << ResidueToDSSPLine(*residue) << std::endl;
    last = residue;
  }
}
