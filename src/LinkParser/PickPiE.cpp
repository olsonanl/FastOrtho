/*
 * PickPiE.cpp
 *
 *  Created on: May 21, 2010
 *      Author: mscott
 */

#include "PickPiE.h"
#include <iostream>
#include <string.h>


//  and island is a continuous block of lines with a common query & subject
// initial read will march through blast file until 1st island is absorbed
// first island will set leftGene & rightGene & eValue
//    liner will absorb first line of next island if one exists
// logic assumes that the e-values from the first island lines
//   in a single query block are never increasing
PickPiE::PickPiE(LineParser *linerSet, GgParser *setGg,
		float setMaxE, double minPercent) {
	// maximum value for 1st eValue in an island
	//     if island starts with a large eValue than
	//     code will skip lines until a new query gene description is
	//     discovered
	maxLogE = setMaxE;

	//  required to support isValid method
	//       computed percent identity value for an island
	//       must be as large as piCut to be considered valid
	piCut = minPercent;

	// allows object to get new lines
	liner = linerSet;

	// allows object to identify query and subject genes in line
	ggParser = setGg;

	//  set try again to true if liner can find first line
	bool tryAgain = liner->getNextLine();

	while (tryAgain) {
		// check to see if block start has known query gene name
		leftGene = ggParser->getGgEntry(liner->getQueryName());
		if (leftGene == NULL) {
			// could not find described query gene so discard line
			tryAgain = liner->getNextLine();
		} else {
			rightGene = ggParser->getGgEntry(liner->getSubjectName());
			if (rightGene == NULL) {
				// could not find described subject gene so discard line
				tryAgain = liner->getNextLine();
			} else {
				// set eValue
				eValue = liner->getValueE();
				if (maxLogE < eValue ) {
					// island is not a candidate
					// advance to a line with  a query name
					//    that does not match leftGene->queryName
					//  set tryAgain if a new query candidate is discovered
					tryAgain = findNextQuery();
				} else {
					// absorb data for computing percent identity
					//    value from first line in island
					percentIdentity = liner->getPercentIdentity();
					matchArea = liner->getMatchArea();
					percentIdentity *= matchArea;

					// finish percent idenity computations and
					//    advance liner to start of next island candidate
					findNextIsland();
					//  first query is a new query
					atNewQuery = true;
					// exit loop
					tryAgain = false;
				}
			}
		}
	}
}

bool PickPiE::isValid() {
	//  make sure percent idenity conditions are met
	return (piCut <= percentIdentity);
}

bool PickPiE::hasBadPi() {
	//  make sure percent idenity conditions are met
	return (percentIdentity < piCut);
}


//  advance line are check to see if query next line
//     has changed
bool PickPiE::setParts() {
	haveParts = liner->getNextLine();
	if (haveParts) {
		lineLeft = liner->getQueryName();
		// check for a change in the query block
		sameLeft =
			(0 == strcmp(leftGene->name, lineLeft));
		lineRight = liner->getSubjectName();
	} else {
		// if all lines are exhausted then
		//   there is no query part to match leftGene
		sameLeft = false;
	}
	return haveParts;
}

// used to get first line of a new query
void PickPiE::startNewQuery() {
	bool tryAgain = false;
	do {
		// make sure left gene name is known
		leftGene = ggParser->getGgEntry(lineLeft);
		if (leftGene == NULL) {
			// move on to next line if query gene is unknown
			tryAgain = liner->getNextLine();
			if (tryAgain) {
				lineLeft = liner->getQueryName();
				lineRight = liner->getSubjectName();
			}
		} else {
			// check to see if right gene is known
			rightGene = ggParser->getGgEntry(lineRight);
			if (rightGene == NULL) {
				// move on to next line if right gene is not known or if
				//   right gene matches left gene
				tryAgain = setParts();
			} else {
				// leftGene & rightGene are known but
				//   e value may be large enough to view island
				//      as a candidate
				eValue = liner->getValueE();
				if (maxLogE < eValue) {
					sameLeft = true;
					// skip remainder of lines with the same query
					tryAgain = findNextQuery();
				} else {
					//  start percent identity calculations for island
					matchArea = liner->getMatchArea();
					percentIdentity =
						matchArea * (liner->getPercentIdentity());
					// advance liner to next island
					//   while completing percent identity calculations
					findNextIsland();
					atNewQuery = true;
					tryAgain = false;
				}
			}
		}
	} while (tryAgain);
}


void PickPiE::findNextIsland() {
	// no more lines will complete next island search
	bool tryAgain = setParts() & sameLeft;
	while (tryAgain) {
		// left genes match so same island is still possible
		tryAgain = (0 == strcmp(rightGene->name, lineRight));
		if (tryAgain) {
			//  add weighted version of percent identity value in
			//    line to percent idenity
			unsigned areaAdd = liner->getMatchArea();
			percentIdentity += areaAdd * (liner->getPercentIdentity());
			//  maintain total area for final
			//    percent identity computation
			matchArea += areaAdd;
			// line was part of island so advance
			tryAgain = setParts() & sameLeft;
		}
	}
	// complete percent identity computations
	percentIdentity /= matchArea;
}

//  liner is pointing at candidate line for a new island
bool PickPiE::haveSameLeft() {
	bool result = false;
	// check to see if next island starts with the same
	//    query as the last
	if (sameLeft) {
		// left gene's match so need to check to see if nextRightGene is
		//    a valid name
		GgEntry* nextRight = ggParser->getGgEntry(lineRight);
		while (nextRight == NULL) {
			// discard line with sameLeft and unknown right query
			setParts();
			// note sameLeft will be set to false if
			//    line data is exhausted
			if (sameLeft) {
				// check to see if this line has a valid subject gene
				nextRight = ggParser->getGgEntry(lineRight);
			} else {
				// left gene name has changed so abandon attempt to
				//   find a valid line with a matching left gene
				nextRight = leftGene;
			}
		}
		if (sameLeft) {
			// found a new island with the same query as the previous
			//    island
			// absorb data from parts
			//    since left gene has not changed
			//    only data concerning right gene needs to be set
			rightGene = nextRight;
			eValue = liner->getValueE();
			if (maxLogE < eValue) {
				// remainder of query block will not produce
				//   candidate islands
				//   so skip to next query block
				findNextQuery();
			} else {
				result = true;
				// begin percent identity computations
				matchArea = liner->getMatchArea();
				percentIdentity = matchArea * (liner->getPercentIdentity());

				//   collect rest of island and finish
				//     percent identity computation
				findNextIsland();
				// indicate that query gene has not changed
				atNewQuery = false;
			}
		}
	}
	return result;
}



bool PickPiE::setNextLine() {
	// use rightGene as flag to indicate that data got extracted
	rightGene = NULL;
	// if all lines are gone exit
	if (haveParts) {
		// if haveSameLeft is true then we got quick extraction from
		//    and next island with the same query
		if ((!haveSameLeft()) && (haveParts)) {
			//  at first line of new query block
			//     still need to get to end of a valid island
			startNewQuery();
			atNewQuery = true;
		}
	}
	return (rightGene != NULL);
}

bool PickPiE::findNextQuery() {
	do {
		setParts();
	} while (sameLeft);
	return haveParts;
}


PickPiE::~PickPiE() {
}
