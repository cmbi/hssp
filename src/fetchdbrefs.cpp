#include "fetchdbrefs.h"

#include "mas.h"
#include "utils.h"

#include <boost/asio.hpp>
#include <boost/foreach.hpp>
#include <boost/format.hpp>
#include <boost/iostreams/device/back_inserter.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/regex.hpp>
#include <zeep/xml/document.hpp>

#include <iostream>

#if defined(_MSC_VER)
#pragma comment ( lib, "libzeep" )
#endif

#define foreach BOOST_FOREACH

using boost::asio::ip::tcp;
using namespace zeep;
namespace io = boost::iostreams;

namespace
{

boost::asio::io_service io_service;
tcp::resolver resolver(io_service);

void FetchHTTPDocument(const std::string& inServer, const std::string& inURL,
                       const std::string& inDb, const std::string& inID,
                       std::string& outHeader, std::string& outDocument)
{
  // This code comes from the http client example code in Boost ASIO

  // Get a list of endpoints corresponding to the server name.
  tcp::resolver::query query(inServer, "http");
  tcp::resolver::iterator endpoint_iterator = resolver.resolve(query);
  tcp::resolver::iterator end;

  // Try each endpoint until we successfully establish a connection.
  tcp::socket socket(io_service);
//  boost::asio::connect(socket, endpoint_iterator);

  boost::system::error_code error = boost::asio::error::host_not_found;
    while (error && endpoint_iterator != end)
    {
    socket.close();
    socket.connect(*endpoint_iterator++, error);
  }
  if (error)
    throw boost::system::system_error(error);

  // Create a soap request in an envelope
  std::string soapRequest = (boost::format(
    "<SOAP-ENV:Envelope xmlns:SOAP-ENV='http://schemas.xmlsoap.org/soap/envelope/'>"
      "<SOAP-ENV:Body>"
        "<ns:GetLinked xmlns:ns='http://mrs.cmbi.umcn.nl/mrsws/search'>"
          "<ns:db>%1%</ns:db>"
          "<ns:id>%2%</ns:id>"
          "<ns:linkedDatabank>pdb</ns:linkedDatabank>"
          "<ns:resultoffset>0</ns:resultoffset>"
          "<ns:maxresultcount>50</ns:maxresultcount>"
        "</ns:GetLinked>"
      "</SOAP-ENV:Body>"
    "</SOAP-ENV:Envelope>") % inDb % inID).str();

  // Form the request. We specify the "Connection: close" header so that the
  // server will close the socket after transmitting the response. This will
  // allow us to treat all data up until the EOF as the content.
  boost::asio::streambuf request;
  std::ostream request_stream(&request);
  request_stream << "POST " << inURL << " HTTP/1.0\r\n";
  request_stream << "Host: " << inServer << "\r\n";
  request_stream << "Accept: */*\r\n";
  request_stream << "User-Agent: mkhssp\r\n";
  request_stream << "Content-Length: " << soapRequest.length() << "\r\n";
  request_stream << "Content-Type: text/xml; charset=\"utf-8\"\r\n";
  request_stream << "Connection: close\r\n\r\n";
  request_stream << soapRequest;

  // Send the request.
  boost::asio::write(socket, request);

  // Read the response status line. The response streambuf will automatically
  // grow to accommodate the entire line. The growth may be limited by passing
  // a maximum size to the streambuf constructor.
  boost::asio::streambuf response;
  boost::asio::read_until(socket, response, "\r\n");

  // Check that response is OK.
  std::istream response_stream(&response);
  std::string http_version;
  response_stream >> http_version;
  unsigned int status_code;
  response_stream >> status_code;
  std::string status_message;
  getline(response_stream, status_message);
  if (!response_stream or http_version.substr(0, 5) != "HTTP/")
    throw mas_exception("Invalid response");

  if (status_code != 200)
    throw mas_exception(
        boost::format("Response returned with status code %d") % status_code);

  // Read the response headers, which are terminated by a blank line.
  boost::asio::read_until(socket, response, "\r\n\r\n");

  // Process the response headers.
  std::string line;

  io::filtering_ostream osh(io::back_inserter(outHeader));

  while (getline(response_stream, line) and line != "\r")
    osh << line << "\n";

  io::filtering_ostream osd(io::back_inserter(outDocument));

  // Write whatever content we already have to output.
  if (response.size() > 0)
    osd << &response;

  // Read until EOF, writing data to output as we go.
//  boost::system::error_code error;
  while (boost::asio::read(socket, response, boost::asio::transfer_at_least(1),
         error))
    osd << &response;

  if (error != boost::asio::error::eof)
    throw boost::system::system_error(error);
}

void FetchHTTPDocument(const std::string& inServer, const std::string& inURL,
                       const std::string& inDb,
                       const std::vector<std::string>& inIDs,
                       std::string& outHeader, std::string& outDocument)
{
  // This code comes from the http client example code in Boost ASIO

  // Get a list of endpoints corresponding to the server name.
  tcp::resolver::query query(inServer, "http");
  tcp::resolver::iterator endpoint_iterator = resolver.resolve(query);
  tcp::resolver::iterator end;

  // Try each endpoint until we successfully establish a connection.
  tcp::socket socket(io_service);
  //  boost::asio::connect(socket, endpoint_iterator);

  boost::system::error_code error = boost::asio::error::host_not_found;
    while (error && endpoint_iterator != end)
    {
    socket.close();
    socket.connect(*endpoint_iterator++, error);
  }
  if (error)
    throw boost::system::system_error(error);

  std::stringstream s;
  s <<
    "<SOAP-ENV:Envelope xmlns:SOAP-ENV='http://schemas.xmlsoap.org/soap/envelope/'>"
      "<SOAP-ENV:Body>"
        "<ns:GetLinkedEx xmlns:ns='http://mrs.cmbi.umcn.nl/mrsws/search'>"
          "<ns:db>" << inDb << "</ns:db>"
          "<ns:linkedDatabank>pdb</ns:linkedDatabank>";

  foreach (std::string id, inIDs)
  s <<       "<ns:id>" << id << "</ns:id>";

  s <<    "</ns:GetLinkedEx>"
      "</SOAP-ENV:Body>"
    "</SOAP-ENV:Envelope>";

  std::string soapRequest = s.str();

  // Form the request. We specify the "Connection: close" header so that the
  // server will close the socket after transmitting the response. This will
  // allow us to treat all data up until the EOF as the content.
  boost::asio::streambuf request;
  std::ostream request_stream(&request);
  request_stream << "POST " << inURL << " HTTP/1.0\r\n";
  request_stream << "Host: " << inServer << "\r\n";
  request_stream << "Accept: */*\r\n";
  request_stream << "User-Agent: mkhssp\r\n";
  request_stream << "Content-Length: " << soapRequest.length() << "\r\n";
  request_stream << "Content-Type: text/xml; charset=\"utf-8\"\r\n";
  request_stream << "Connection: close\r\n\r\n";
  request_stream << soapRequest;

  // Send the request.
  boost::asio::write(socket, request);

  // Read the response status line. The response streambuf will automatically
  // grow to accommodate the entire line. The growth may be limited by passing
  // a maximum size to the streambuf constructor.
  boost::asio::streambuf response;
  boost::asio::read_until(socket, response, "\r\n");

  // Check that response is OK.
  std::istream response_stream(&response);
  std::string http_version;
  response_stream >> http_version;
  unsigned int status_code;
  response_stream >> status_code;
  std::string status_message;
  getline(response_stream, status_message);
  if (!response_stream or http_version.substr(0, 5) != "HTTP/")
    throw mas_exception("Invalid response");

  if (status_code != 200)
    throw mas_exception(
        boost::format("Response returned with status code %d") % status_code);

  // Read the response headers, which are terminated by a blank line.
  boost::asio::read_until(socket, response, "\r\n\r\n");

  // Process the response headers.
  std::string line;

  io::filtering_ostream osh(io::back_inserter(outHeader));

  while (getline(response_stream, line) and line != "\r")
    osh << line << "\n";

  io::filtering_ostream osd(io::back_inserter(outDocument));

  // Write whatever content we already have to output.
  if (response.size() > 0)
    osd << &response;

  // Read until EOF, writing data to output as we go.
//  boost::system::error_code error;
  while (boost::asio::read(socket, response, boost::asio::transfer_at_least(1),
                           error))
    osd << &response;

  if (error != boost::asio::error::eof)
    throw boost::system::system_error(error);
}

}

void FetchPDBReferences(const std::string& inBaseURL, const std::string& inDb,
                        const std::string& inID,
                        std::vector<std::string>& outReferences)
{
  static const boost::regex re("http://([^/]+)(/.*)?");

  boost::smatch m;
  if (not boost::regex_match(inBaseURL, m, re))
    throw mas_exception("Invalid base url for FetchPDBReferences");

  outReferences.clear();

  try
  {
    std::string httpHeader, httpDocument;
    FetchHTTPDocument(m[1], inBaseURL, inDb, inID, httpHeader, httpDocument);

    xml::document doc(httpDocument);

    auto hits = doc.find("//hits/id");
    foreach (auto hit, hits)
      outReferences.push_back(hit->content());
  }
  catch (...)
  {

  }
}

void FetchPDBReferences(const std::string& inBaseURL, const std::string& inDb,
  std::map<std::string,std::vector<std::string>>& ioReferences)
{
  static const boost::regex re("http://([^/]+)(/.*)?");

  boost::smatch m;
  if (not boost::regex_match(inBaseURL, m, re))
    throw mas_exception("Invalid base url for FetchPDBReferences");

  try
  {
    std::vector<std::string> ids;
    foreach (auto r, ioReferences)
      ids.push_back(r.first);

    std::string httpHeader, httpDocument;
    FetchHTTPDocument(m[1], inBaseURL, inDb, ids, httpHeader, httpDocument);

    xml::document doc(httpDocument);

    foreach (auto r, doc.find("//response"))
    {
      auto id = r->find_first("id");
      if (not id)
        continue;

      std::vector<std::string> links;
      foreach (auto link, r->find("linked"))
        links.push_back(link->content());

      ioReferences[id->content()] = links;
    }
  }
  catch (...)
  {

  }
}
