/*
 * PickE.cpp
 *
 *  Created on: May 19, 2010
 *      Author: mscott
 */

#include "PickE.h"
#include <iostream>
#include <string.h>


//  and island is a continuous block of lines with a common query & subject
// initial read will march through blast file until 1st island is absorbed
// first island will set leftGene & rightGene & eValue
//    liner will absorb first line of next island if one exists
// logic assumes that the e-values from the first island lines
//   in a single query block are never increasing
PickE::PickE(LineParser *linerSet, GgParser *setGg, float setMaxE) {
	lineLeft = NULL;
	lineRight = NULL;
	//  allows skipping islands that with large e-values
	maxLogE = setMaxE;
	//  used to collect links from m8 or bpo type file
	//   this will filter links where the query and subject are identical
	liner = linerSet;
	// required to converting gene names into indices
	ggParser = setGg;
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
					// entire block is invalid
					// advance to a line with  a query name
					//    that does not match leftGene->queryName
					tryAgain = findNextQuery();
				} else {
					// found an island with valid gene names and e-value
					// process remainder of island
					findNextIsland();
					// first query is a new query
					atNewQuery = true;
					// exit loop
					tryAgain = false;
				}
			}
		}
	}
}


bool PickE::setParts() {
	// use liner to get parts
	haveParts = liner->getNextLine();
	if (haveParts) {
		lineLeft = liner->getQueryName();
		// check for a change in the query block
		sameLeft = (0 == strcmp(leftGene->name, lineLeft));
		lineRight = liner->getSubjectName();
	} else {
		// if all lines are exhausted then
		//   there is no query part to match leftGene
		sameLeft = false;
	}
	return haveParts;
}

// In PickE large e-value is only failure method but
//    but these islands are never returned
bool PickE::isValid() {
	return true;
}

// In PickE large e-value is only failure method but
//    but these islands are never returned
bool PickE::hasBadPi() {
	return false;
}

// at first line in a new query block
void PickE::startNewQuery() {
	bool tryAgain = false;
	do {
		// make sure left gene name is known
		leftGene = ggParser->getGgEntry(lineLeft);
		if (leftGene == NULL) {
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
				//   e value may be invalid
				eValue = liner->getValueE();
				if (maxLogE < eValue) {
					sameLeft = true;
					tryAgain = findNextQuery();
				} else {
					// advance liner to next island
					findNextIsland();
					atNewQuery = true;
					tryAgain = false;
				}
			}
		}
	} while (tryAgain);
}

void PickE::findNextIsland() {
	// if end of lines then end of island
	bool tryAgain = setParts() & sameLeft;;
	while (tryAgain) {
		// left genes match so same island is still possible
		tryAgain = (0 == strcmp(rightGene->name, lineRight));
		if (tryAgain) {
			// line was part of island so advance
				tryAgain = setParts() & sameLeft;
		}
	}
}

// liner contains a line contents
//    sameLeft indicates whether query in line matches leftGene
//    nextRight
bool PickE::haveSameLeft() {
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
			// note: sameLeft will be set to false if
			//    line data is exhausted
			if (sameLeft) {
				// check to see if this line has a valid subject gene
				nextRight = ggParser->getGgEntry(lineRight);
			} else {
				// left gene name has changed so abandon attempt to
				//   find a valid line with a matching left gene
				//  by placing a non-NULL value in nextRight
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
				// remainder of query block is invalid
				//   so skip to next query block
				findNextQuery();
			} else {
				result = true;
				// indicate that query gene has not changed
				atNewQuery = false;
				//   collect rest of island
				findNextIsland();
			}
		}
	}
	return result;
}

bool PickE::setNextLine() {
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
		}
	}
	return (rightGene != NULL);

}

bool PickE::findNextQuery() {
	do {
		setParts();
	} while (sameLeft);
	return haveParts;
}


PickE::~PickE() {
}
