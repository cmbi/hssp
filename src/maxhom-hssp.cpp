//  Copyright Maarten L. Hekkelman, Radboud University 2008.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include "MRS.h"

#include <wait.h>

#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/filesystem/convenience.hpp>
#include <boost/filesystem/fstream.hpp>
#include <boost/foreach.hpp>
#define foreach BOOST_FOREACH

#include "mas.h"
#include "blast.h"
#include "dssp.h"
#include "structure.h"
#include "maxhom-hssp.h"

using namespace std;
namespace ba = boost::algorithm;

namespace maxhom
{

void WriteToFD(
	int						inFD,
	string					inText)
{
	inText += '\n';
	
	for (;;)
	{
		int r = write(inFD, inText.c_str(), inText.length());
		if (r == -1 and errno == EAGAIN)
			continue;
		if (r != static_cast<int>(inText.length()))
			THROW(("Failed to write maxhom command"));
		break;
	}		 
}

void GetHSSPForHitsAndDSSP(
	CDatabankPtr			inDatabank,
	const string&			inMaxHom,
	const string&			inPDBID,
	const vector<uint32>&	inHits,
	const string&			inDSSP,
	int						inMaxAlign,
	ostream&				outHSSP)
{
	int threshold = 5;
	
	HUuid uuid;
	
	fs::path rundir("/tmp/hssp-2/");
	rundir /= boost::lexical_cast<string>(uuid);
	fs::create_directories(rundir);
	
	// create 'hits' file containing swiss-prot entries for all hits
	fs::ofstream hits(rundir / "hits");
	for (vector<uint32>::const_iterator hit = inHits.begin(); hit != inHits.end(); ++hit)
		hits << inDatabank->GetDocument(*hit);
	hits.close();
	
	// write out the dssp file
	fs::ofstream dssp(rundir / (inPDBID + ".dssp"));
	dssp << inDSSP;
	dssp.close();
	
	// start a maxhom
	int ifd[2];
	
	if (pipe(ifd) != 0)
		THROW(("Failed to create pipe: %s", strerror(errno)));
	
	int pid = fork();
	
	if (pid == -1)
	{
		close(ifd[0]);
		close(ifd[1]);
		
		THROW(("fork failed: %s", strerror(errno)));
	}
	
	if (pid == 0)	// the child process (will be maxhom)
	{
		fs::current_path(rundir);
		
		setpgid(0, 0);
		
		dup2(ifd[0], STDIN_FILENO);
		close(ifd[0]);
		close(ifd[1]);
		
		int fd = open("maxhom.log", O_CREAT | O_RDWR, 0666);
		dup2(fd, STDOUT_FILENO);
		dup2(fd, STDERR_FILENO);
		
		char* argv[] = {
			const_cast<char*>(inMaxHom.c_str()),
			nil
		};
		
		(void)execve(inMaxHom.c_str(), argv, environ);
	}
	
	close(ifd[0]);
	int fd = ifd[1];
	
	// OK, so now we pass maxhom the commands and wait for it to exit
	WriteToFD(fd, "COMMAND NO");
	WriteToFD(fd, "BATCH");
	WriteToFD(fd, (boost::format("PID: %d") % getpid()).str());
	WriteToFD(fd, "SEQ_1 " + inPDBID + ".dssp");
	WriteToFD(fd, "SEQ_2 hits");
	WriteToFD(fd, "2_PROFILES NO");
	WriteToFD(fd, "METRIC LACHLAN");
	WriteToFD(fd, "NORM_PROFILE DISABLED");
	WriteToFD(fd, "MEAN_PROFILE ignored");
	WriteToFD(fd, "FACTOR_GAPS ignored");
	WriteToFD(fd, "SMIN -0.50");
	WriteToFD(fd, "SMAX 1.00");
	WriteToFD(fd, "GAP_OPEN 3.00");
	WriteToFD(fd, "GAP_ELONG 0.10");
	WriteToFD(fd, "WEIGHT1 YES");
	WriteToFD(fd, "WEIGHT2 NO");
	WriteToFD(fd, "WAY3-ALIGN NO");
	WriteToFD(fd, "INDEL_1 YES");
	WriteToFD(fd, "INDEL_2 YES");
	WriteToFD(fd, "RELIABILITY NO");
	WriteToFD(fd, "FILTER_RANGE 10.0");
	WriteToFD(fd, "NBEST 1");
	WriteToFD(fd, (boost::format("MAXALIGN %d") % inMaxAlign).str());
	WriteToFD(fd, (boost::format("THRESHOLD FORMULA%+d") % threshold).str());
	WriteToFD(fd, "SORT DISTANCE");
	WriteToFD(fd, "HSSP out.hssp");
	WriteToFD(fd, "SAME_SEQ_SHOW YES");
	WriteToFD(fd, "SUPERPOS NO");
	WriteToFD(fd, "PDB_PATH OPTION DISABLED");
	WriteToFD(fd, "PROFILE_OUT NO");
	WriteToFD(fd, "STRIP_OUT NO");
	WriteToFD(fd, "LONG_OUT NO");
	WriteToFD(fd, "DOT_PLOT NO");
	WriteToFD(fd, "RUN");
	
	int status;
	waitpid(pid, &status, 0);
	
	if (status != 0)
		THROW(("maxhom exited with status %d", status));
	
	// OK, got it! Read in the result and exit
	fs::ifstream result(rundir / "out.hssp");

	// update the SEQBASE line while copying over the data
	string dbVersion = inDatabank->GetVersion();
	
	for (;;)
	{
		string line;
		getline(result, line);
		if (result.eof())
			break;
		
		if (ba::starts_with(line, "SEQBASE    "))
		{
			line.erase(11, string::npos);
			line += dbVersion;
		}
		
		outHSSP << line << endl;
	}
}

void CreateHSSP(
	CDatabankPtr		inDatabank,
	const string&		inMaxHom,
	MProtein&			inProtein,
	ostream&			outHSSP)
{
	stringstream dssp;
	WriteDSSP(inProtein, dssp);

	vector<string> seqs;
	inProtein.GetSequences(back_inserter(seqs));
	seqs.erase(unique(seqs.begin(), seqs.end()), seqs.end());
	
	vector<uint32> hits;
	foreach (const string& seq, seqs)
	{
		vector<uint32> h1;
		BlastSequence(inDatabank, seq, h1);
		sort(h1.begin(), h1.end());
		
		vector<uint32> h2;
		set_union(h1.begin(), h1.end(), hits.begin(), hits.end(), back_inserter(h2));
		swap(h2, hits);
	}
	
	if (hits.empty())
		THROW(("No blast hits found for %s", inProtein.GetID().c_str()));
	
	string hssp;
	GetHSSPForHitsAndDSSP(inDatabank, inMaxHom, inProtein.GetID(), hits, dssp.str(), 2000, outHSSP);
}

}
