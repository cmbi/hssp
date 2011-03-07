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

	auto_ptr<CDocIterator> data(new CDbAllDocIterator(inDatabank.get()));
	CBlast blast(inSequence, matrix, wordsize, expect, filter, gapped, gapOpen, gapExtend, maxhits);
	
	if (blast.Find(*inDatabank, *data, threads))
	{
		CBlastHitList hits(blast.Hits());
		
		foreach (const CBlastHit& hit, hits)
			outHits.push_back(hit.DocumentNr());
	}
}

