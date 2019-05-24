// Copyright Maarten L. Hekkelman, Radboud University 2008.
// Copyright Coos Baakman, Jon Black, Wouter G. Touw & Gert Vriend, Radboud university medical center 2015.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "fasta.h"
#include "hssp-nt.h"
#include "mas.h"
#include "structure.h"
#include "utils.h"

#include <boost/algorithm/string.hpp>
#include <boost/filesystem/fstream.hpp>
#include <boost/foreach.hpp>
#include <boost/iostreams/copy.hpp>
#include <boost/iostreams/device/back_inserter.hpp>

#ifdef HAVE_LIBBZ2
#include <boost/iostreams/filter/bzip2.hpp>
#endif
#ifdef HAVE_LIBZ
#include <boost/iostreams/filter/gzip.hpp>
#endif

#include <boost/iostreams/filter/zlib.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/program_options.hpp>
#include <boost/program_options/config.hpp>

#if defined(_MSC_VER)
#include <conio.h>
#endif

namespace ba = boost::algorithm;
namespace fs = boost::filesystem;
namespace io = boost::iostreams;
namespace po = boost::program_options;

#define foreach BOOST_FOREACH
// int VERBOSE = 0;


int main(int argc, char* argv[])
{
#if P_UNIX
  // enable the dumping of cores to enable postmortem debugging
  rlimit l;
  if (getrlimit(RLIMIT_CORE, &l) == 0)
  {
    l.rlim_cur = l.rlim_max;
    if (l.rlim_cur == 0 or setrlimit(RLIMIT_CORE, &l) < 0)
      std::cerr << "Failed to set rlimit" << std::endl;
  }
#endif

  fs::path outfilename;

  try
  {
    po::options_description desc("mkhssp options");
    desc.add_options()
      ("help,h", "Display help message")
      ("input,i", po::value<std::string>(), "PDB ID or input PDB file (.pdb), mmCIF file (.cif/.mcif), or fasta file (.fa/.fasta), optionally compressed by gzip (.gz) or bzip2 (.bz2)")
      ("output,o", po::value<std::string>(), "Output file, optionally compressed by gzip (.gz) or bzip2 (.bz2). Use 'stdout' to output to screen")
      ("databank,d", po::value<std::vector<std::string>>(), "Databank to use (can be specified multiple times)")
      ("threads,a",  po::value<uint32>(), "Number of threads (default is maximum)")
//      ("use-seqres", po::value<bool>(), "Use SEQRES chain instead of chain based on ATOM records (values are true of false, default is true)")
      ("min-length", po::value<uint32>(), "Minimal chain length (default = 25)")
      ("fragment-cutoff", po::value<float>(), "Minimal alignment length as fraction of chain length (default = 0.75)")
      ("gap-open,O", po::value<float>(), "Gap opening penalty (default is 30.0)")
      ("gap-extend,E", po::value<float>(), "Gap extension penalty (default is 2.0)")
      ("threshold", po::value<float>(), "Homology threshold adjustment (default = 0.05)")
      ("max-hits,m", po::value<uint32>(), "Maximum number of hits to include (default = 5000)")
#ifdef HAVE_LIBZEEP
      ("fetch-dbrefs", "Fetch DBREF records for each UniProt ID")
#endif
      ("verbose,v", "Verbose mode")
      ("version", "Show version number and citation info");

    po::positional_options_description p;
    p.add("input", 1);
    p.add("output", 2);

    po::variables_map vm;
    po::store(po::command_line_parser(argc, argv).options(desc).positional(p).run(), vm);

    try
    {
      fs::path home = get_home();
      if (fs::exists(home / ".mkhssprc"))
      {
        fs::ifstream rc(home / ".mkhssprc");
        po::store(po::parse_config_file(rc, desc), vm);
      }
    }
    catch (const mas_exception& e)
    {
      // TODO: Log a debug message when a proper logging library is used.
    }

    po::notify(vm);

    if (vm.count("version")>0)
    {
      std::cout << "mkhssp version " << PACKAGE_VERSION << std::endl
         << std::endl
         << "If you use HSSP, please cite: " << std::endl
         << "Touw WG, Baakman C, Black J, te Beek TA, Krieger E, Joosten RP & Vriend G." << std::endl
         << "A series of PDB-related databanks for everyday needs." << std::endl
         << "Nucleic Acids Res. (2015) 43, D364-D368. doi: 10.1093/nar/gku1028." << std::endl
         << std::endl
         << "The original HSSP reference is: " << std::endl
         << "Sander C & Schneider R." << std::endl
         << "Database of homology-derived protein structures and the structural meaning of sequence alignment." << std::endl
         << "Proteins (1991) 9, 56-68. doi: 10.1002/prot.340090107." << std::endl;
      exit(0);
    }

    if (vm.count("help") or not vm.count("input") or vm.count("databank") == 0)
    {
      std::cerr << desc << std::endl;
      exit(1);
    }

    VERBOSE = vm.count("verbose");

    std::vector<fs::path> databanks;
    std::vector<std::string> dbs = vm["databank"].as<std::vector<std::string>>();
    foreach (std::string db, dbs)
    {
      databanks.push_back(db);
      if (not fs::exists(databanks.back()))
        throw mas_exception(boost::format("Databank %s does not exist") % db);
    }

//    bool useSeqRes = true;
//    if (vm.count("use-seqres"))
//      useSeqRes = vm["use-seqres"].as<bool>();

    uint32 minlength = 25;
    if (vm.count("min-length"))
      minlength= vm["min-length"].as<uint32>();

    uint32 maxhits = 5000;
    if (vm.count("max-hits"))
      maxhits= vm["max-hits"].as<uint32>();

    float gapOpen = 30;
    if (vm.count("gap-open"))
      gapOpen = vm["gap-open"].as<float>();

    float gapExtend = 2;
    if (vm.count("gap-extend"))
      gapExtend = vm["gap-extend"].as<float>();

    float threshold = HSSP::kThreshold;
    if (vm.count("threshold"))
      threshold = vm["threshold"].as<float>();

    float fragmentCutOff = HSSP::kFragmentCutOff;
    if (vm.count("fragment-cutoff"))
      fragmentCutOff = vm["fragment-cutoff"].as<float>();

#ifdef HAVE_LIBZEEP
    bool fetchDbRefs = vm.count("fetch-dbrefs") > 0;
#else
    bool fetchDbRefs = false;
#endif

    uint32 threads = boost::thread::hardware_concurrency();
    if (vm.count("threads"))
      threads = vm["threads"].as<uint32>();
    if (threads < 1)
      threads = 1;

    // what input to use
    std::string input = vm["input"].as<std::string>();
    io::filtering_stream<io::input> in;
    std::ifstream infile(input.c_str(), std::ios_base::in | std::ios_base::binary);
    if (not infile.is_open())
      throw std::runtime_error("Error opening input file");

#ifdef HAVE_LIBBZ2
    if (ba::ends_with(input, ".bz2"))
    {
      in.push(io::bzip2_decompressor());
      input.erase(input.length() - 4, std::string::npos);
    }
#endif
#ifdef HAVE_LIBZ
    if (ba::ends_with(input, ".gz"))
    {
      in.push(io::gzip_decompressor());
      input.erase(input.length() - 3, std::string::npos);
    }
#endif
    in.push(infile);

    // Where to write our HSSP file to:
    // either to cout or an (optionally compressed) file.
    std::ofstream outfile;
    io::filtering_stream<io::output> out;

    if (vm.count("output") and vm["output"].as<std::string>() != "stdout")
    {
      outfilename = fs::path(vm["output"].as<std::string>());
      outfile.open(outfilename.c_str(), std::ios_base::out|std::ios_base::trunc|std::ios_base::binary);

      if (not outfile.is_open())
        throw std::runtime_error("could not create output file");

#ifdef HAVE_LIBBZ2
      if (ba::ends_with(outfilename.string(), ".bz2"))
        out.push(io::bzip2_compressor());
#endif
#ifdef HAVE_LIBZ
      if (ba::ends_with(outfilename.string(), ".gz"))
        out.push(io::gzip_compressor());
#endif

      out.push(outfile);
    }
    else
      out.push(std::cout);

    // if input file is a FastA file, we process it differently
    if (ba::ends_with(input, ".fa") or ba::ends_with(input, ".fasta"))
    {
      std::vector<MProtein*> proteins = read_proteins_from_fasta(in);
      for (auto& p: proteins)
      {
        try
        {
          HSSP::CreateHSSP(*p, databanks, maxhits, minlength, gapOpen,
                           gapExtend, threshold, fragmentCutOff, threads,
                           fetchDbRefs, out);
        }
        catch (const std::exception& e)
        {
          std::cerr << "Creating HSSP for " << p->GetID() << " failed: "
                    << e.what() << std::endl;
        }

        delete p;
      }
    }
    else
    {
      // read protein and calculate the secondary structure
      MProtein a;

      if (ba::ends_with(input, ".cif") or ba::ends_with(input, ".mcif"))
        a.ReadmmCIF(in);
      else
        a.ReadPDB(in);

      a.CalculateSecondaryStructure();

      // create the HSSP file
      HSSP::CreateHSSP(a, databanks, maxhits, minlength,
        gapOpen, gapExtend, threshold, fragmentCutOff, threads, fetchDbRefs, out);
    }
  }
  catch (const std::exception& e)
  {
    std::cerr << e.what() << std::endl;

    try
    {
      if (not outfilename.empty() and fs::exists(outfilename))
        fs::remove(outfilename);
    }
    catch (...) {}

    exit(1);
  }

//#if defined(_MSC_VER) && ! NDEBUG
//  std::cerr << "Press any key to quit application ";
//  char ch = _getch();
//  std::cerr << std::endl;
//#endif

  return 0;
}

