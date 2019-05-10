#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <boost/algorithm/string.hpp>
#include <boost/filesystem/convenience.hpp>
#include <boost/filesystem/fstream.hpp>
#include <boost/foreach.hpp>
#include <boost/iostreams/copy.hpp>

#ifdef HAVE_LIBBZ2
#include <boost/iostreams/filter/bzip2.hpp>
#endif
#ifdef HAVE_LIBZ
#include <boost/iostreams/filter/gzip.hpp>
#endif

#include <boost/format.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/program_options.hpp>
#include <boost/program_options/config.hpp>
#include <boost/ptr_container/ptr_vector.hpp>
#include <boost/regex.hpp>

#include "mas.h"
#include "utils.h"

#include "hssp-convert-3to1.h"

namespace ba = boost::algorithm;
namespace io = boost::iostreams;
namespace fs = boost::filesystem;
namespace po = boost::program_options;

int VERBOSE = 0;

int main(int argc, char* const argv[])
{
  try
  {
    po::options_description desc("hsspconv options");
    desc.add_options()
      ("help,h", "Display help message")
      ("input,i", po::value<std::string>(), "Input PDB file (or PDB ID)")
      ("output,o",
       po::value<std::string>(),
       "Output file, use 'stdout' to output to screen")
      ("version", "Show version number");

    po::positional_options_description p;
    p.add("input", 1);
    p.add("output", 2);

    po::variables_map vm;
    po::store(po::command_line_parser(
          argc, argv).options(desc).positional(p).run(), vm);

    //fs::path home = get_home();
    //if (fs::exists(home / ".mkhssprc"))
    //{
    //  fs::ifstream rc(home / ".mkhssprc");
    //  po::store(po::parse_config_file(rc, desc), vm);
    //}

    po::notify(vm);

    if(vm.count("version")>0)
    {
      std::cout << "hssp converter version "<< PACKAGE_VERSION << std::endl;
      exit(0);
    }

    if (vm.count("help"))
    {
      std::cerr << desc << std::endl;
      exit(1);
    }
    io::filtering_stream<io::input> in;
    fs::ifstream ifs;

    if (vm.count("input") == 0)
      in.push(std::cin);
    else
    {
      fs::path input = vm["input"].as<std::string>();
      ifs.open(input, std::ios::binary);

      if (not ifs.is_open())
        throw mas_exception(
            boost::format("Could not open input file '%s'") % input);

#ifdef HAVE_LIBBZ2
      if (input.extension() == ".bz2")
        in.push(io::bzip2_decompressor());
#endif
#ifdef HAVE_LIBZ
      if (input.extension() == ".gz")
        in.push(io::gzip_decompressor());
#endif
      in.push(ifs);
    }

    if (vm.count("output") == 0)
      ConvertHsspFile(in, std::cout);
    else
    {
      fs::path output = vm["output"].as<std::string>();

      fs::ofstream ofs(output, std::ios::binary);
      if (not ofs.is_open())
        throw mas_exception(
            boost::format("Could not open output file '%s'") % output);

      io::filtering_stream<io::output> out;

#ifdef HAVE_LIBBZ2
      if (output.extension() == ".bz2")
        out.push(io::bzip2_compressor());
#endif
#ifdef HAVE_LIBZ
      if (output.extension() == ".gz")
        out.push(io::gzip_compressor());
#endif
      out.push(ofs);

      ConvertHsspFile(in, out);
    }
  }
  catch (const std::exception& e)
  {
    std::cerr << e.what() << std::endl;
    exit(1);
  }
  catch (...)
  {
    std::cerr << "Unknown exception" << std::endl;
    exit(1);
  }

  return 0;
}
