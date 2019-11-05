#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Fasta


#include "fasta.h"

#include <boost/test/unit_test.hpp>

#include <sstream>


BOOST_AUTO_TEST_SUITE(test_mkhssp_suite)

BOOST_AUTO_TEST_CASE(test_read_proteins_from_fasta_single)
{
  std::stringstream ss;
  ss << ">test" << std::endl;
  ss << "TTCCPSIVARSNFNVCRLPGTPEAICATYTGCIIIPGATCPGDYAN" << std::endl;
  auto proteins = read_proteins_from_fasta(ss);

  BOOST_CHECK_EQUAL(proteins.size(), unsigned(1));
}

BOOST_AUTO_TEST_CASE(test_read_proteins_from_fasta_multiple)
{
  std::stringstream ss;
  ss << ">test1" << std::endl;
  ss << "TTCCPSIVARSNFNVCRLPGTPEAICATYTGCIIIPGATCPGDYAN" << std::endl;
  ss << ">test2" << std::endl;
  ss << "TTCCPSIVARSNFNVCRLPGTPEAICATYTGCIIIPGATCPGDYAN" << std::endl;
  ss << ">test3" << std::endl;
  ss << "TTCCPSIVARSNFNVCRLPGTPEAICATYTGCIIIPGATCPGDYAN" << std::endl;
  auto proteins = read_proteins_from_fasta(ss);

  BOOST_CHECK_EQUAL(proteins.size(), unsigned(3));
}

BOOST_AUTO_TEST_SUITE_END()
