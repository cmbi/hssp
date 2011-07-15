//  Copyright Maarten L. Hekkelman, Radboud University 2008.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include "MRS.h"

#include <boost/foreach.hpp>
#define foreach BOOST_FOREACH
#include <boost/thread.hpp>

#include "CBlast.h"

#include "mas.h"
#include "blast.h"
#include "structure.h"

using namespace std;

int BLAST_THREADS = 0;

void BlastSequence(
	CDatabankPtr		inDatabank,
	const string&		inSequence,
	vector<uint32>&		outHits)
{
	float expect = 1.0;
	bool filter = true, gapped = true;
	int wordsize = 3, gapOpen = 11, gapExtend = 1, threads, maxhits = 1500;
	string matrix = "BLOSUM62";
	
	threads = BLAST_THREADS;
	if (threads < 1)
		threads = boost::thread::hardware_concurrency();

//	auto_ptr<CDocIterator> data(new CDbAllDocIterator(inDatabank.get()));
//	CBlast blast(inSequence, matrix, wordsize, expect, filter, gapped, gapOpen, gapExtend, maxhits);

	CBlastResult* results = inDatabank->PerformBlastSearch(
		inSequence, "blastp", matrix, wordsize, expect, filter, gapped, gapOpen, gapExtend, maxhits);
	
	if (results != nil)
	{
		foreach (const CBlastHit& hit, results->hits)
			outHits.push_back(hit.documentNr);
		
		delete results;
	}
}

void BlastProtein(
	CDatabankPtr		inDatabank,
	const MProtein&		inProtein,
	vector<uint32>&		outHits)
{
	foreach (const MChain* chain, inProtein.GetChains())
	{
		char chainID = chain->GetChainID();
		entry e(1, inProtein.GetID() + chainID);
		inProtein.GetSequence(chainID, e);

		BlastSequence(inDatabank, decode(e.m_seq), outHits);
	}
}
