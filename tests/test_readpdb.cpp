#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE ReadPDB

#include "structure.h"

#include <boost/test/unit_test.hpp>

#include <sstream>

namespace xssp { namespace test { } }

using namespace xssp::test;

using std::istream;
using std::istringstream;
using std::string;

namespace xssp {
  namespace test {

    class readpdb_fixture {
    protected:
      ~readpdb_fixture() noexcept = default;

      /// \brief Read an istream of PDB data into an MProtein and return it
      ///
      /// This may become redundant if MProtein's interface is changed to provide this directly
      static MProtein ReadPDBIntoMProtein(std::istream &argIs ///< The istream of PDB data to parse
                                          ) {
        MProtein theProtein;
        theProtein.ReadPDB( argIs );
        return theProtein;
      }

      /// \brief Get the number of residues in the specified MProtein
      ///
      /// This may become redundant if the MProtein's interface is changed to provide this directly
      static uint32 getNrOfResidues(const MProtein &argProtein ///< The MProtein to query
                                    ) {
        uint32 nrOfResidues, nrOfChains, nrOfSSBridges, nrOfIntraChainSSBridges, nrOfHBonds;
        uint32 nrOfHBondsPerDistance[11] = {};
        argProtein.GetStatistics(nrOfResidues, nrOfChains, nrOfSSBridges, nrOfIntraChainSSBridges, nrOfHBonds, nrOfHBondsPerDistance);
        return nrOfResidues;
      }

      /// \brief Perform Boost Test assertions that the specified residue is present in the specified protein
      ///
      /// This checks that getting the residue doesn't throw, which is the current way of indicating absence.
      /// It also checks the result isn't a nullptr, in case that behvaiour is changed.
      static void checkResidueIsPresent(const MProtein &argProtein,            ///< The protein to query
                                        const string   &argChainID,            ///< The chain ID of the residue to check
                                        const uint16   &argSeqNumber,          ///< The seq number of the residue to check
                                        const string   &argInsertionCode = " " ///< The insertion code of the residue to check
                                                                               ///<   (with  adefault of " " for empty insert code)
                                        ) {
        BOOST_REQUIRE_NO_THROW( argProtein.GetResidue( argChainID, argSeqNumber, argInsertionCode )            );
        BOOST_CHECK           ( argProtein.GetResidue( argChainID, argSeqNumber, argInsertionCode ) != nullptr );
      }
    };

  }
}



BOOST_FIXTURE_TEST_SUITE(readpdb_test_suite, readpdb_fixture)



BOOST_AUTO_TEST_CASE(parses_four_normal_residues)
{
  // Given: this raw data
  istringstream raw_pdb_data_ss{ R"(ATOM      1  N   PRO A   1       3.069   2.269 -37.532  1.00 45.60           N  
ATOM      2  CA  PRO A   1       3.988   2.650 -36.447  1.00 44.31           C  
ATOM      3  C   PRO A   1       3.746   1.828 -35.194  1.00 38.99           C  
ATOM      4  O   PRO A   1       2.578   1.409 -35.069  1.00 44.26           O  
ATOM      8  N   PRO A   2       4.679   1.563 -34.283  1.00 33.16           N  
ATOM      9  CA  PRO A   2       4.222   0.910 -33.032  1.00 28.69           C  
ATOM     10  C   PRO A   2       3.183   1.806 -32.357  1.00 24.26           C  
ATOM     11  O   PRO A   2       3.215   3.030 -32.488  1.00 26.34           O  
ATOM     15  N   GLY A   3       2.280   1.100 -31.665  1.00 20.53           N  
ATOM     16  CA  GLY A   3       1.219   1.833 -30.986  1.00 17.75           C  
ATOM     17  C   GLY A   3       1.770   2.570 -29.801  1.00 15.03           C  
ATOM     18  O   GLY A   3       2.952   2.505 -29.471  1.00 17.11           O  
ATOM     19  N   PRO A   4       0.904   3.267 -29.124  1.00 14.96           N  
ATOM     20  CA  PRO A   4       1.397   4.003 -27.950  1.00 14.56           C  
ATOM     21  C   PRO A   4       1.644   3.066 -26.798  1.00 12.09           C  
ATOM     22  O   PRO A   4       1.184   1.937 -26.714  1.00 12.02           O  
)" };

  // When: parsing the data into a MProtein
  const MProtein theProtein = ReadPDBIntoMProtein( raw_pdb_data_ss );

  // Then: the resulting MProtein should have 4 residues: A:1, A:2, A:3, A:4
  BOOST_CHECK_EQUAL( getNrOfResidues( theProtein ), 4 );
  checkResidueIsPresent( theProtein, "A", 1 );
  checkResidueIsPresent( theProtein, "A", 2 );
  checkResidueIsPresent( theProtein, "A", 3 );
  checkResidueIsPresent( theProtein, "A", 4 );
}



// Test the issue documented in https://github.com/cmbi/xssp/issues/79
//
// Residue 124's ATOM records were being ignored because the code didn't
// adequately reset after rejecting the residue 123's ATOM records.
BOOST_AUTO_TEST_CASE(includes_valid_altlocn_residue_after_rejecting_altlocn_residue)
{
  // Given: this raw data
  istringstream raw_pdb_data_ss{ R"(ATOM   1032  N   MET A 120      23.127   5.241  -5.234  1.00  9.79           N  
ATOM   1033  CA  MET A 120      24.061   4.468  -6.042  1.00  9.61           C  
ATOM   1034  C   MET A 120      24.811   3.451  -5.187  1.00  9.38           C  
ATOM   1035  O   MET A 120      25.054   2.306  -5.620  1.00 10.06           O  
ATOM   1065  N  BGLN A 123      22.567   0.603  -4.409  1.00 11.11           N  
ATOM   1066  CA BGLN A 123      22.431  -0.200  -5.625  1.00 11.15           C  
ATOM   1067  C  BGLN A 123      23.706  -0.950  -5.973  1.00 10.67           C  
ATOM   1068  O  BGLN A 123      23.706  -1.700  -6.955  1.00 11.70           O  
ATOM   1074  N  ALYS A 124      24.772  -0.666  -5.315  0.64  9.24           N  
ATOM   1075  N  BLYS A 124      24.774  -0.724  -5.226  0.36 10.75           N  
ATOM   1076  CA ALYS A 124      26.058  -1.322  -5.576  0.64  9.02           C  
ATOM   1077  CA BLYS A 124      26.116  -1.237  -5.480  0.36 10.92           C  
ATOM   1078  C  ALYS A 124      26.612  -0.964  -6.962  0.64  9.77           C  
ATOM   1079  C  BLYS A 124      26.599  -0.944  -6.904  0.36  9.69           C  
ATOM   1080  O  ALYS A 124      27.304  -1.768  -7.599  0.64 11.47           O  
ATOM   1081  O  BLYS A 124      27.229  -1.776  -7.565  0.36 11.86           O  
ATOM   1092  N   ARG A 125      26.324   0.274  -7.353  1.00  9.58           N  
ATOM   1093  CA  ARG A 125      26.849   0.839  -8.587  1.00 10.21           C  
ATOM   1094  C   ARG A 125      28.105   1.626  -8.216  1.00  9.67           C  
ATOM   1095  O   ARG A 125      28.075   2.852  -8.074  1.00 10.39           O  
)"
  };

  // When: parsing the data into a MProtein
  const MProtein theProtein = ReadPDBIntoMProtein( raw_pdb_data_ss );

  // Then: the resulting MProtein should have at least 3 residues
  // (but accept more in case A:123 is accepted in the future)
  // including A:120, A:124 and A:125
  BOOST_CHECK_GE( getNrOfResidues( theProtein ), 3 );
  checkResidueIsPresent( theProtein, "A", 120 );
  checkResidueIsPresent( theProtein, "A", 124 );
  checkResidueIsPresent( theProtein, "A", 125 );
}



BOOST_AUTO_TEST_SUITE_END()
