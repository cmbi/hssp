#include "mas.h"

#include <boost/format.hpp>
#include <boost/regex.hpp>
#include <boost/asio.hpp>
#include <boost/iostreams/device/back_inserter.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/foreach.hpp>
#define foreach BOOST_FOREACH

#include <zeep/xml/document.hpp>

#include "fetchdbrefs.h"
#include "utils.h"

#if defined(_MSC_VER)
#pragma comment ( lib, "libzeep" )
#endif

using namespace std;
using boost::asio::ip::tcp;
using namespace zeep;
namespace io = boost::iostreams;

namespace
{

boost::asio::io_service io_service;
tcp::resolver resolver(io_service);

void FetchHTTPDocument(const string& inServer, const string& inURL,
	const string& inDb, const string& inID, string& outHeader, string& outDocument)
{
	// This code comes from the http client example code in Boost ASIO

	// Get a list of endpoints corresponding to the server name.
	tcp::resolver::query query(inServer, "http");
	tcp::resolver::iterator endpoint_iterator = resolver.resolve(query);
	tcp::resolver::iterator end;

	// Try each endpoint until we successfully establish a connection.
	tcp::socket socket(io_service);
//	boost::asio::connect(socket, endpoint_iterator);

	boost::system::error_code error = boost::asio::error::host_not_found;
    while (error && endpoint_iterator != end)
    {
		socket.close();
		socket.connect(*endpoint_iterator++, error);
	}
	if (error)
		throw boost::system::system_error(error);

	// Create a soap request in an envelope
	string soapRequest = (boost::format(
		"<SOAP-ENV:Envelope xmlns:SOAP-ENV='http://schemas.xmlsoap.org/soap/envelope/'>"
			"<SOAP-ENV:Body>"
				"<ns:GetLinked xmlns:ns='http://mrs.cmbi.ru.nl/mrsws/search'>"
					"<ns:db>%1%</ns:db>"
					"<ns:id>%2%</ns:id>"
					"<ns:linkedDatabank>pdb</ns:linkedDatabank>"
				"</ns:GetLinked>"
			"</SOAP-ENV:Body>"
		"</SOAP-ENV:Envelope>") % inDb % inID).str();

	// Form the request. We specify the "Connection: close" header so that the
	// server will close the socket after transmitting the response. This will
	// allow us to treat all data up until the EOF as the content.
	boost::asio::streambuf request;
	ostream request_stream(&request);
	request_stream << "POST " << inURL << " HTTP/1.0\r\n";
	request_stream << "Host: " << inServer << "\r\n";
	request_stream << "Accept: */*\r\n";
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
	istream response_stream(&response);
	string http_version;
	response_stream >> http_version;
	unsigned int status_code;
	response_stream >> status_code;
	string status_message;
	getline(response_stream, status_message);
	if (response_stream == false or http_version.substr(0, 5) != "HTTP/")
		throw mas_exception("Invalid response");

	if (status_code != 200)
		throw mas_exception(boost::format("Response returned with status code %d") % status_code);

	// Read the response headers, which are terminated by a blank line.
	boost::asio::read_until(socket, response, "\r\n\r\n");

	// Process the response headers.
	string line;

	io::filtering_ostream osh(io::back_inserter(outHeader));

	while (getline(response_stream, line) and line != "\r")
		osh << line << "\n";

	io::filtering_ostream osd(io::back_inserter(outDocument));

	// Write whatever content we already have to output.
	if (response.size() > 0)
		osd << &response;

	// Read until EOF, writing data to output as we go.
//	boost::system::error_code error;
	while (boost::asio::read(socket, response, boost::asio::transfer_at_least(1), error))
		osd << &response;

	if (error != boost::asio::error::eof)
		throw boost::system::system_error(error);
}

}

void FetchPDBReferences(const string& inBaseURL, const string& inDb, const string& inID,
	vector<string>& outReferences)
{
	static const boost::regex re("http://([^/]+)(/.*)?");

	outReferences.clear();

	try
	{
		boost::smatch m;
		if (not boost::regex_match(inBaseURL, m, re))
			throw mas_exception("Invalid base url for FetchPDBReferences");
		
		string httpHeader, httpDocument;
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
