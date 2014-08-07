// Copyright Maarten L. Hekkelman, Radboud University 2008-2011.
//   Distributed under the Boost Software License, Version 1.0.
//       (See accompanying file LICENSE_1_0.txt or copy at
//             http://www.boost.org/LICENSE_1_0.txt)

#include "blast.h"
#include "dssp.h"
#include "hssp-nt.h"
#include "mas.h"
#include "structure.h"
#include "utils.h"
#include "version.h"

#include <boost/algorithm/string.hpp>
#include <boost/filesystem/convenience.hpp>
#include <boost/filesystem/fstream.hpp>
#include <boost/filesystem/operations.hpp>
#include <boost/foreach.hpp>
#include <boost/format.hpp>
#include <boost/iostreams/copy.hpp>
#include <boost/iostreams/device/back_inserter.hpp>
#include <boost/iostreams/filter/newline.hpp>
#include <boost/iostreams/filter/zlib.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/program_options.hpp>
#include <boost/regex.hpp>
#include "zeep/config.hpp"
#include "zeep/server.hpp"

#include <pwd.h>
#include <signal.h>
#include <sys/resource.h>


#define foreach BOOST_FOREACH
#define HSSPSOAP_PID_FILE "/var/run/hsspsoap.pid"
#define HSSPSOAP_LOG_FILE "/var/log/hsspsoap.log"

namespace ba = boost::algorithm;
namespace fs = boost::filesystem;
namespace io = boost::iostreams;
namespace po = boost::program_options;

// embedded resources
extern char _binary_rsrc_index_html_start[];
extern char _binary_rsrc_index_html_size[];
extern char _binary_rsrc_error_html_start[];
extern char _binary_rsrc_error_html_size[];

// globals
//int VERBOSE;
uint32 gNrOfThreads;

void GetDSSPForSequence(const std::string& inSequence, std::string& outDSSP)
{
  io::filtering_ostream out(io::back_inserter(outDSSP));

  out << "==== Secondary Structure Definition by the program DSSP, updated CMBI version by ElmK / April 1,2000 ==== DATE=28-MAY-2010     ." << std::endl
    << "REFERENCE" << std::endl
    << "HEADER                                                        9UNK" << std::endl
    << "COMPND" << std::endl
    << "SOURCE" << std::endl
    << "AUTHOR" << std::endl
    << boost::format("%5.5d  1  0  0  0") % inSequence.length() << std::endl
    << "  #" << std::endl;

  // And now fill in the rest
  int n = 1;
  std::string::const_iterator aa;
  for (aa = inSequence.begin(); aa != inSequence.end(); ++aa, ++n)
    out << boost::format("%5.d%5.d A %c") % n % n % char(toupper(*aa))
        << std::endl;
}

void GetPDBFileFromPayload(const std::string& payload, std::string& pdb,
                           fs::path& file)
{
  // make streams
  io::filtering_istream in;
  in.push(io::newline_filter(io::newline::posix));
  in.push(boost::make_iterator_range(payload));

  // get the boundary
  std::string boundary;
  getline(in, boundary);

  // parse fields until we've got the data 'pdb'
  std::string name;

  for (;;)
  {
    // we just read a boundary, what follows are header fields
    std::string line;
    getline(in, line);

    // sanity check
    if (line.empty() and in.eof())
      throw mas_exception("Unexpected end of file");

    // skip header fields until we have "Content-Disposition: form-data"
    if (ba::starts_with(line, "Content-Disposition: form-data"))
    {
      static const boost::regex nre("\\bname=\\\"([^\"]+)\\\"");
      boost::smatch m;

      if (boost::regex_search(line, m, nre))
        name = m[1];
      else
        name = "undef";  // really??

      static const boost::regex fre("\\bfilename=\\\"([^\"]+)\\\"");
      if (boost::regex_search(line, m, fre))
        file = m[1].str();

      continue;
    }

    if (not line.empty())  // any other field header
      continue;

    // the data, read until we hit the next boundary

    std::string data;
    io::filtering_ostream out(io::back_inserter(pdb));

    for (;;)
    {
      getline(in, line);

      if (line.empty() and in.eof())
        throw mas_exception("Unexpected end of file");

      if (ba::starts_with(line, boundary))
        break;

      if (name == "pdb" or name == "pdbfile")
        out << line << std::endl;
    }

    // check to see if we're done
    if ((name == "pdb" or name == "pdbfile") and pdb.length() > 2)
      break;

    if (line.substr(boundary.length(), 2) == "--")
      break;
  }
}

// the data types used in our communication with the outside world
// are wrapped in a namespace.

class hssp_server : public zeep::server
{
  public:
    hssp_server(const std::vector<fs::path>& inDatabank);

  virtual void handle_request(const zeep::http::request& req,
                              zeep::http::reply& rep);

  virtual void GetDSSPForPDBFile(const std::string& pdbfile,
                                 std::string& dssp);
  virtual void GetHSSPForPDBFile(const std::string& pdbfile,
                                 std::string& hssp);
  virtual void GetHSSPForSequence(const std::string& sequence,
                                  std::string& hssp);

  std::vector<fs::path>
          mDatabank;
};

hssp_server::hssp_server(const std::vector<fs::path>& inDatabank)
  : zeep::server("http://www.cmbi.ru.nl/hsspsoap", "hsspsoap"),
    mDatabank(inDatabank)
{
  const char* kGetDSSPForPDBFileParameterNames[] = {"pdbfile", "dssp"};
  register_action("GetDSSPForPDBFile", this, &hssp_server::GetDSSPForPDBFile,
                  kGetDSSPForPDBFileParameterNames);

  const char* kGetHSSPForPDBFileParameterNames[] = {"pdbfile", "hssp"};
  register_action("GetHSSPForPDBFile", this, &hssp_server::GetHSSPForPDBFile,
                  kGetHSSPForPDBFileParameterNames);

  const char* kGetHSSPForSequenceParameterNames[] = {"sequence", "hssp"};
  register_action("GetHSSPForSequence", this, &hssp_server::GetHSSPForSequence,
                  kGetHSSPForSequenceParameterNames);
}

void hssp_server::handle_request(const zeep::http::request& req,
                                 zeep::http::reply& rep)
{
  bool handled = false;

  std::string uri = req.uri;

  // strip off the http part including hostname and such
  if (ba::starts_with(uri, "http://"))
  {
    std::string::size_type s = uri.find_first_of('/', 7);
    if (s != std::string::npos)
      uri.erase(0, s);
  }

  // now make the path relative to the root
  while (uri.length() > 0 and uri[0] == '/')
    uri.erase(uri.begin());

  try
  {
    if (req.method == "GET" and
        (uri.empty() or ba::starts_with(uri, "index.htm")))
    {
      rep.set_content(std::string(_binary_rsrc_index_html_start,
                                  _binary_rsrc_index_html_size), "text/html");
      handled = true;
    }
    else if (req.method == "POST")
    {
      if (ba::starts_with(uri, "PDB2DSSP") or
          ba::starts_with(uri, "PDB2HSSP") )
      {
        std::string pdb;
        fs::path file;

        GetPDBFileFromPayload(req.payload, pdb, file);

        if (file.empty() and pdb.length() > 66)
          file = pdb.substr(62, 4) + ".pdb";

        std::string result;
        if (ba::starts_with(uri, "PDB2DSSP"))
        {
          GetDSSPForPDBFile(pdb, result);
          file.replace_extension(".dssp");
        }
        else
        {
          GetHSSPForPDBFile(pdb, result);
          file.replace_extension(".hssp");
        }

        rep.set_content(result, "text/plain");
        rep.set_header("Content-disposition",
          (boost::format("attachement; filename=\"%1%\"") % file).str());

        handled = true;
      }
      else if (ba::starts_with(uri, "SEQ2HSSP") )
      {
        std::string::size_type p = req.payload.find("seq=");
        if (p == std::string::npos)
          throw mas_exception("Missing sequence parameters");

        std::string seq = req.payload.substr(p + 4);
        seq = zeep::http::decode_url(seq);

        std::string result;
        GetHSSPForSequence(seq, result);

        rep.set_content(result, "text/plain");
        rep.set_header("Content-disposition",
          (boost::format("attachement; filename=\"%1%\"") % "hssp-for-sequence").str());

        handled = true;
      }
    }
  }
  catch (const std::exception& e)
  {
    //mrsrc::rsrc rsrc("error.html");
    //std::string error(rsrc.data(), rsrc.size());

    std::string error(_binary_rsrc_error_html_start, _binary_rsrc_error_html_size);

    ba::replace_first(error, "#ERRSTR", e.what());

    rep.set_content(error, "text/html");
    handled = true;
  }

  if (not handled)
    zeep::server::handle_request(req, rep);
}

void hssp_server::GetDSSPForPDBFile(const std::string& pdbfile,
                                    std::string& dssp)
{
  // create a protein
  io::filtering_istream in(boost::make_iterator_range(pdbfile));
  MProtein a;
  a.ReadPDB(in);

  // then calculate the secondary structure
  a.CalculateSecondaryStructure();

  io::filtering_ostream out(io::back_inserter(dssp));
  WriteDSSP(a, out);
}

void hssp_server::GetHSSPForPDBFile(const std::string& pdbfile,
                                    std::string& hssp)
{
  io::filtering_istream in;
  in.push(io::newline_filter(io::newline::posix));
  in.push(boost::make_iterator_range(pdbfile));

  // OK, we've got the file, now create a protein
  MProtein a;
  a.ReadPDB(in);

  // then calculate the secondary structure
  a.CalculateSecondaryStructure();

  // finally, create the HSSP
  io::filtering_ostream out(io::back_inserter(hssp));
  HSSP::CreateHSSP(a, mDatabank, 5000, 25, 30, 2,
    HSSP::kThreshold, HSSP::kFragmentCutOff,
    gNrOfThreads, false, out);
}

void hssp_server::GetHSSPForSequence(const std::string& sequence,
                                     std::string& hssp)
{
  io::filtering_ostream out(io::back_inserter(hssp));
  HSSP::CreateHSSP(sequence, mDatabank, 5000, 25, 30, 2, HSSP::kThreshold,
                   HSSP::kFragmentCutOff, gNrOfThreads, false, out);
}

// --------------------------------------------------------------------
//
//  Daemonize
//

void Daemonize(const std::string& inUser)
{
  int pid = fork();

  if (pid == -1)
  {
    std::cerr << "Fork failed" << std::endl;
    exit(1);
  }

  if (pid != 0)
    _exit(0);

  if (setsid() < 0)
  {
    std::cerr << "Failed to create process group: " << strerror(errno)
              << std::endl;
    exit(1);
  }

  // it is dubious if this is needed:
  signal(SIGHUP, SIG_IGN);

  // fork again, to avoid being able to attach to a terminal device
  pid = fork();

  if (pid == -1)
    std::cerr << "Fork failed" << std::endl;

  if (pid != 0)
    _exit(0);

  // write our pid to the pid file
  std::ofstream pidFile(HSSPSOAP_PID_FILE);
  pidFile << getpid() << std::endl;
  pidFile.close();

  if (chdir("/") != 0)
  {
    std::cerr << "Cannot chdir to /: " << strerror(errno) << std::endl;
    exit(1);
  }

  if (inUser.length() > 0)
  {
    struct passwd* pw = getpwnam(inUser.c_str());
    if (pw == NULL or setuid(pw->pw_uid) < 0)
    {
      std::cerr << "Failed to set uid to " << inUser << ": "
                << strerror(errno) << std::endl;
      exit(1);
    }
  }

  // close stdin
  close(STDIN_FILENO);
  open("/dev/null", O_RDONLY);
}

// --------------------------------------------------------------------
//
//  OpenLogFile
//

void OpenLogFile()
{
  // open the log file
  int fd = open(HSSPSOAP_LOG_FILE, O_CREAT|O_APPEND|O_RDWR, 0644);
  if (fd < 0)
  {
    std::cerr << "Opening log file " HSSPSOAP_LOG_FILE " failed" << std::endl;
    exit(1);
  }

  // redirect stdout and stderr to the log file
  dup2(fd, STDOUT_FILENO);
  dup2(fd, STDERR_FILENO);
}

// --------------------------------------------------------------------
//
//  main
//

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

  po::options_description desc("Options");
  desc.add_options()
    ("help,h", "Display help message")
    ("address,a", po::value<std::string>(),  "address to bind to")
    ("port,p", po::value<uint16>(),  "port to bind to")
    ("location,l", po::value<std::string>(),  "location advertised in wsdl")
    ("user,u", po::value<std::string>(),  "user to run as")
    ("threads,a", po::value<uint32>(),  "number of threads to use")
    ("databank", po::value<std::vector<std::string>>(),
     "Databank(s) to use (default = /data/fasta/uniprot_sprot.fasta and /data/fasta/uniprot_trembl.fasta)")
    ("no-daemon,D", "do not fork a daemon")
    ("verbose,v", "Verbose mode")
    ("version", "Show version number")
    ;

  std::string location = "http://www.cmbi.ru.nl/hsspsoap/wsdl";
  std::string address = "0.0.0.0";
  std::string user = "nobody";

  uint16 port = 10334;

  bool daemon = true;

  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, desc), vm);
  po::notify(vm);

  if (vm.count("version") > 0)
  {
    std::cout << "hsspsoap version " << XSSP_VERSION << std::endl;
    exit(0);
  }

  if (vm.count("help") or vm.count("databank") == 0)
  {
    std::cout << desc << std::endl;
    exit(1);
  }

  if (vm.count("address"))
    address = vm["address"].as<std::string>();

  if (vm.count("location"))
    location = vm["location"].as<std::string>();

  if (vm.count("port"))
    port = vm["port"].as<uint16>();

  if (vm.count("user"))
    user = vm["user"].as<std::string>();

  gNrOfThreads = boost::thread::hardware_concurrency();
  if (vm.count("threads"))
    gNrOfThreads = vm["threads"].as<uint32>();
  if (gNrOfThreads < 1)
    gNrOfThreads = 1;

  if (vm.count("verbose"))
    VERBOSE = 1;

  std::vector<fs::path> databanks;
  std::vector<std::string> dbs = vm["databank"].as<std::vector<std::string>>();
  foreach (std::string db, dbs)
  {
    databanks.push_back(db);
    if (not fs::exists(databanks.back()))
      throw mas_exception(boost::format("Databank %s does not exist") % db);
  }

  if (vm.count("no-daemon"))
    daemon = false;

  if (daemon)
  {
    OpenLogFile();
    Daemonize(user);
  }

#ifndef _MSC_VER
    sigset_t new_mask, old_mask;
    sigfillset(&new_mask);
    pthread_sigmask(SIG_BLOCK, &new_mask, &old_mask);
#endif

  // create server
  hssp_server server(databanks);
  server.bind(address, port);

  if (not location.empty())
    server.set_location(location);

    boost::thread_group t;
    t.create_thread(boost::bind(&hssp_server::run, &server, 4));

#ifndef _MSC_VER
    pthread_sigmask(SIG_SETMASK, &old_mask, 0);

  // Wait for signal indicating time to shut down.
  sigset_t wait_mask;
  sigemptyset(&wait_mask);
  sigaddset(&wait_mask, SIGINT);
  sigaddset(&wait_mask, SIGQUIT);
  sigaddset(&wait_mask, SIGTERM);
  pthread_sigmask(SIG_BLOCK, &wait_mask, 0);
  int sig = 0;
  sigwait(&wait_mask, &sig);

  server.stop();
#endif

  t.join_all();

  return 0;
}
