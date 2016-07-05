/*
 * BpoAllParser.cpp
 *
 *  Created on: May 21, 2010
 *      Author: mscott
 */

#include "BpoAllParser.h"
#include <string.h>
#include <iostream>


BpoAllParser::BpoAllParser(GgParser* ggSet, const char *bpoFile, char splitSet,
							double maxWeight, unsigned *indexSettings) {
	splitter = splitSet;
	// make sure minColumnCount is not too large
	minColumnCount = 0;
	// set column locations from inputs while updatine minColumnCount
	setIndex(QUERY_AT, indexSettings, &queryAt);
	setIndex(QUERY_LENGTH, indexSettings, &queryLength);
	setIndex(SUBJECT_AT, indexSettings, &subjectAt);
	setIndex(SUBJECT_LENGTH, indexSettings, &subjectLength);
	setIndex(E_VALUE_AT, indexSettings, &eValueAt);
	setIndex(IDENTITY_AT, indexSettings, &identityAt);
	setIndex(MAP_COLUMN, indexSettings, &mapColumn);

	// minColumnCount is now set at max index
	//  increment since column indexes are 0 based
	++minColumnCount;

	//  set value to use for log10(0.0)
	minLogE = -maxWeight;

	//  prepare to march through the file
	puller = new LinePull(0x10000, bpoFile);

	// save ggParser so gene lengths can be set from bpo lines
	ggParser = ggSet;
	if (ggParser != NULL) {
		// prepare ggParser to accept gene lengths
		ggParser->createLengthsRoom();
	}
	// set appropriate intialization values
	nextMap = NULL;
	lastQueryName = NULL;
	lastSubjectName = NULL;
}

void BpoAllParser::setIndex(unsigned findAt, unsigned *settings,
							unsigned *toSet) {
	*toSet = settings[findAt];
	if (minColumnCount < *toSet) {
		minColumnCount = *toSet;
	}
}

char * BpoAllParser::getQueryName() {
	return (parts[queryAt]);
}

char *BpoAllParser::getSubjectName() {
	return (parts[subjectAt]);
}

float BpoAllParser::getValueE() {
	return (eReturn);
}

char* BpoAllParser::getTextE() {
	return (parts[eValueAt]);
}

unsigned BpoAllParser::getQueryStart() {
	return (queryStart);
}

unsigned BpoAllParser::getSubjectStart() {
	return (subjectStart);
}

unsigned BpoAllParser::getQueryEnd() {
	return (queryEnd);
}

unsigned BpoAllParser::getSubjectEnd() {
	return (subjectEnd);
}

unsigned BpoAllParser::getMatchArea() {
	//  match area is not provided by bpo lines but
	//     the bpo creation process should block multi-line islands
	return (1);
}

float BpoAllParser::getPercentIdentity() {
	return (atof(parts[identityAt]));
}

char *setInt(char *start, char end, unsigned *value) {
	char *result = NULL;
	// see if any text follows last termination position
	if ('\0' != *start) {
		// search for character that terminates next int value
		result = strchr(start, end);
		if (result != NULL) {
			// start should now point a integer text
			//   terminate this text so its value can be extracted with
			//   atoi()
			*result = '\0';
			++result;
			*value = atoi(start);
		}
	}
	return result;
}

bool BpoAllParser::getMapValues(char *mapText) {
	//  assume nothing is available
	nextMap = NULL;
	if (mapText != NULL) {
		// skip past index of map part
		mapText = strchr(mapText, ':');
		if (mapText != NULL) {
			++mapText;
			//  next part in query start with - termination character
			mapText = setInt(mapText, '-', &queryStart);
			if (mapText != NULL) {
				//  next part is query end with : termination character
				mapText = setInt(mapText, ':', &queryEnd);
				if (mapText != NULL) {
					// next part is subject start with - termination character
					mapText = setInt(mapText, '-', &subjectStart);
					if (mapText != NULL) {
						//  last part is subject end with . termination
						nextMap = setInt(mapText, '.', &subjectEnd);
					}
				}
			}
		}
	}
	return (nextMap != NULL);
}


bool BpoAllParser::getNextLine() {
	// check to see if the line should produce multi-line returns due
	//    to a multi-part mapping column
	if (!getMapValues(nextMap)) {
		// get all of the parts in the next line from the bpo file
		unsigned partCount;
		do {
			parts = puller->getCharSplitParts(splitter, partCount);
			if (parts == NULL) {
				// bpo file is exhausted so set set partCount to exit loop
				partCount = minColumnCount;
			} else if (minColumnCount <= partCount) {
				// all required parts are present
				// skip line if query == subject or if
				//   invalid data in mapping column
				if ((0 == strcmp(parts[queryAt], parts[subjectAt]))
								||
						(!getMapValues(parts[mapColumn]))){
					partCount = 0;
				} else {
					// decode e-value column
					eReturn = getValueLog(minLogE, parts[eValueAt]);
					if (ggParser != NULL){
						// save gene length data in ggParser
						lastQueryName =
							ggParser->addBpoLength(lastQueryName,
												   parts[queryAt],
												   parts[queryLength]);
						lastSubjectName =
							ggParser->addBpoLength(lastSubjectName,
												   parts[subjectAt],
												   parts[subjectLength]);
					}
				}
			}
		} while (partCount < minColumnCount);
	}
	return (parts != NULL);
}

BpoAllParser::~BpoAllParser() {
}
