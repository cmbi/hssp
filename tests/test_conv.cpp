#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Convert

#include <boost/test/unit_test.hpp>
#include <sstream>
#include <fstream>

#include "hssp-convert-3to1.h"
#include "utils.h"

BOOST_AUTO_TEST_SUITE(test_hsspconv_suite)

BOOST_AUTO_TEST_CASE(test_conv_4s1h)
{
    std::ifstream in;
    in.open("test-data/wrong1.sto");

    std::ostringstream out;

    BOOST_CHECK_THROW(
        ConvertHsspFile(in, out),
        mas_exception
    );

    in.close();
}

BOOST_AUTO_TEST_SUITE_END()
