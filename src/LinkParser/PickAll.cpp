/*
 * PickAll.cpp
 *
 *  Created on: May 21, 2010
 *      Author: mscott
 */

#include "PickAll.h"
#include <string.h>

#define START_AT 0
#define END_AT 1
#define MAP_PARTS 2

PickAll::PickAll(LineParser *linerSet, GgParser *setGg,
		float setMaxE, double minPercentIdentity, double minMatchCover) {
	// maximum value for 1st eValue in an island
	//     if island starts with a large eValue than
	//     code will skip lines until a new query gene description is
	//     discovered
	maxLogE = setMaxE;

	//  required to support isValid method
	//       computed percent identity value for an island
	//       must be as large as piCut to be considered valid
	piCut = minPercentIdentity;

	// also required for computing percent match
	matchCut = minMatchCover;

	// required for percent matching computations
	geneLengths = setGg->geneLengths;

	liner = linerSet;
	ggParser = setGg;

	// used to support percent matching computations
	maps = new Blocks(MAP_PARTS * sizeof(unsigned));
	bool tryAgain = liner->getNextLine();
	while (tryAgain) {
		leftGene = ggParser->getGgEntry(liner->getQueryName());
		if (leftGene == NULL) {
			tryAgain = liner->getNextLine();
		} else {
			rightGene = ggParser->getGgEntry(liner->getSubjectName());
			if (rightGene == NULL) {
				tryAgain = liner->getNextLine();
			} else {
				// set eValue and validity flag
				eValue = liner->getValueE();
				if (maxLogE < eValue) {
					// block is not a candidate
					// advance to a line with  a query name
					//    that does not match leftGene->queryName
					//  set tryAgain if a new query candidate is discovered
					tryAgain = findNextQuery();

				} else {
					// set for percent match computations
					matchArea = liner->getMatchArea();
					// collect and store area used in 1st matching
					unsigned *mapStore = (unsigned *)(maps->getAdd());
					// load queryLength and subjectLength
					//   to determine minimum
					queryLength = geneLengths[leftGene->order];
					subjectLength = geneLengths[rightGene->order];
					if (queryLength <= subjectLength) {
						// track mapping in query area since
						//    it is not larger
						setMapSpan(liner->getQueryStart(),
								   liner->getQueryEnd(), mapStore);
						smallerLength = queryLength;
					} else {
						// tracking mapping in subject area since
						//    query area is larger
						setMapSpan(liner->getSubjectStart(),
								   liner->getSubjectEnd(), mapStore);
						smallerLength = subjectLength;
					}
					// begin percent identity computations
					percentIdentity =
						matchArea * (liner->getPercentIdentity());
					// finish percent identity and match computations and
					//    advance liner to start of next island candidate
					findNextIsland();
					//  first query is a new query
					atNewQuery = true;
					// exit loop
					tryAgain = false;
				}
			}
		}
		sameLeft = false;
	}
}

void PickAll::setMapSpan(unsigned start, unsigned end, unsigned *store) {
	if (start < end) {
		store[START_AT] = start;
		store[END_AT] = end + 1;
	} else {
		store[START_AT] = end;
		store[END_AT] = start + 1;
	}
}

bool PickAll::isValid() {
	//  check percent identity and matching coverage conditions
	return ((piCut <= percentIdentity) && (matchCut <= percentMatch));
}

bool PickAll::hasBadPi() {
	//  check percent identity and matching coverage conditions
	return (percentIdentity < piCut);
}


bool PickAll::setParts() {
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
void PickAll::startNewQuery() {
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
				// set eValue and validity flag
				eValue = liner->getValueE();
				if (maxLogE < eValue) {
					// block is not a candidate
					// advance to a line with  a query name
					//    that does not match leftGene->queryName
					//  set tryAgain if a new query candidate is discovered
					tryAgain = findNextQuery();
				} else {
					// collect and store area used in 1st matching
					matchArea = liner->getMatchArea();
					unsigned *mapStore = (unsigned *)(maps->getAdd());
					// load queryLength and subjectLength
					//   to determine minimum
					queryLength = geneLengths[leftGene->order];
					subjectLength = geneLengths[rightGene->order];
					if (queryLength <= subjectLength) {
						// track mapping in query area since
						//    it is not larger
						setMapSpan(liner->getQueryStart(),
								   liner->getQueryEnd(), mapStore);
						smallerLength = queryLength;
					} else {
						// tracking mapping in subject area since
						//    query area is larger
						setMapSpan(liner->getSubjectStart(),
								   liner->getSubjectEnd(), mapStore);
						smallerLength = subjectLength;
					}

					// begin percent idenity computations
					percentIdentity =
						matchArea * (liner->getPercentIdentity());
					// finish percent idenity computations and
					//    advance lner to start of next island candidate
					findNextIsland();
					//  first query is a new query
					atNewQuery = true;
					// exit loop
					tryAgain = false;
				}
			}
		}
	} while (tryAgain);
}

// sort that simplifies computation of area involved in matching
int compareMatching(const char *left, const char *right) {
	unsigned *leftInts = (unsigned *)left;
	unsigned *rightInts = (unsigned *)right;
	int result = 0;
	if (leftInts[START_AT] < rightInts[START_AT]) {
		--result;
	} else if (rightInts[START_AT] < leftInts[START_AT]) {
		++result;
	} else {
		if (leftInts[END_AT] < rightInts[END_AT]) {
			--result;
		} else if (rightInts[END_AT] < leftInts[END_AT]) {
			++result;
		}
	}
	return result;
}


void PickAll::findNextIsland() {
	// if end of lines then end of island
	bool tryAgain = setParts() & sameLeft;
	while (tryAgain) {
		// left genes match so extra line status is still possible
		tryAgain = (0 == strcmp(rightGene->name, lineRight));
		if (tryAgain) {
			// subject gene also matches
			unsigned areaAdd = liner->getMatchArea();
			percentIdentity += areaAdd * (liner->getPercentIdentity());
			matchArea += areaAdd;
			unsigned *mapStore = (unsigned *)(maps->getAdd());
			if (queryLength <= subjectLength) {
				// track mapping in query area since
				//    it is not larger
				setMapSpan(liner->getQueryStart(),
						liner->getQueryEnd(), mapStore);
			} else {
				// tracking mapping in subject area since
				//    query area is larger
				setMapSpan(liner->getSubjectStart(),
						   liner->getSubjectEnd(), mapStore);
			}
			// both genes match so advance to next line
			tryAgain = setParts() & sameLeft;
		}
	}
	percentIdentity /= matchArea;

	// sort collected mapping sections to compute
	//    area of gene involved in mapping
	maps->sort(compareMatching);
	BlockMarch* marcher = maps->getFullMarch();
	unsigned *reader = (unsigned *)(marcher->getEntry());
	unsigned matchTop = 0;
	unsigned cover = 0;

	while (reader != NULL) {
		unsigned place = reader[START_AT];
		if (matchTop <= place) {
			// next map area is above collected map area
			// so entire interval contributes to cover
			cover += reader[END_AT] - reader[START_AT];
			//  new top is top of next interval
			matchTop = reader[END_AT];
		} else {
			// find top of next interval
			place = reader[END_AT];
			// no more work if top is not extended
			if (matchTop < place) {
				// extend top and add extension amount to cover
				cover += place - matchTop;
				matchTop = place;
			}
		}
		reader = (unsigned *)(marcher->getNextEntry());
	}
	delete(marcher);
	// clear maps for next island
	maps->empty();
	// use values to test match coverage
	percentMatch = cover;
	percentMatch /= smallerLength;
}


bool PickAll::haveSameLeft() {
	bool result = false;
	// check to see if next island starts with the same
	//    query as the last
	if (sameLeft) {
		// left gene's match so need to check to see if nextRightGene is
		//    a valid name
		GgEntry* nextRight = ggParser->getGgEntry(lineRight);
		while (nextRight == NULL) {
			setParts();
			// note sameLeft will be set to false if
			//    line data is exhausted
			if (sameLeft) {
				// have another chance of a new island with the
				//    same query as the previous island
				nextRight = ggParser->getGgEntry(lineRight);
			} else {
				// left gene name has changed so abandon attempt to
				//   find a valid line with a matching left gene
				nextRight = leftGene;
			}
		}
		if (sameLeft) {
			// found a new island with the same query as the previous island
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
				matchArea = liner->getMatchArea();

				subjectLength = geneLengths[rightGene->order];
				unsigned *mapStore = (unsigned *)(maps->getAdd());
				if (queryLength <= subjectLength) {
					// track mapping in query area since
					//    it is not larger
					setMapSpan(liner->getQueryStart(),
							   liner->getQueryEnd(), mapStore);
					smallerLength = queryLength;
				} else {
					// tracking mapping in subject area since
					//    query area is larger
					setMapSpan(liner->getSubjectStart(),
							   liner->getSubjectEnd(), mapStore);
					smallerLength = subjectLength;
				}
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

bool PickAll::setNextLine() {
	// use rightGene as flag to indicate that data got extracted
	rightGene = NULL;
	if (haveParts) {
		// try to extract data from line with matching left gene
		if ((!haveSameLeft()) && (haveParts)) {
			// failed to find island with old query
			//   but have a line with a new query
			startNewQuery();
			atNewQuery = true;
		}
	}
	return (rightGene != NULL);
}

bool PickAll::findNextQuery() {
	do {
		setParts();
	} while (sameLeft);
	return haveParts;
}


PickAll::~PickAll() {
	delete(maps);
}
