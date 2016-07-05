/*
 * CharSplitParser.cpp
 *
 *  Created on: Apr 28, 2010
 *      Author: mscott
 */

#include "CharSplitParser.h"
#include <iostream>
#include <math.h>
#include <stdlib.h>
#include <string.h>

CharSplitParser::CharSplitParser(char splitSet, const char *tabFile,
		float maxWeight, unsigned *indexSettings) {
	splitter = splitSet;
	puller = new LinePull(0x10000, tabFile);
	// provide value for log10(0.0)
	minLogE = -maxWeight;

	// insure a valid initial value for minSplitCount
	minSplitCount = 0;

	// set column locations and track maximum column position in minSplitCount
	setIndex(QUERY_AT, indexSettings, &queryAt);
	setIndex(QUERY_START, indexSettings, &queryStartAt);
	setIndex(QUERY_END, indexSettings, &queryEndAt);

	setIndex(SUBJECT_AT, indexSettings, &subjectAt);
	setIndex(SUBJECT_START, indexSettings, &subjectStartAt);
	setIndex(SUBJECT_END, indexSettings, &subjectEndAt);

	setIndex(E_VALUE_AT, indexSettings, &eValueAt);
	setIndex(IDENTITY_AT, indexSettings, &identityAt);
	setIndex(LENGTH_AT, indexSettings, &areaAt);

	// adjust minSplitCount to conform with 0 based column indexing
	++minSplitCount;
}

void CharSplitParser::setIndex(unsigned findAt, unsigned *settings,
							   unsigned *toSet) {
	*toSet = settings[findAt];
	if (minSplitCount < *toSet) {
		minSplitCount = *toSet;
	}
}

char * CharSplitParser::getQueryName() {
	return (parts[queryAt]);
}

char *CharSplitParser::getSubjectName() {
	return (parts[subjectAt]);
}

float CharSplitParser::getValueE() {
	return (getValueLog(minLogE, parts[eValueAt]));
}

char* CharSplitParser::getTextE() {
	return (parts[eValueAt]);
}

unsigned CharSplitParser::getQueryStart() {
	return (atoi(parts[queryStartAt]));
}

unsigned CharSplitParser::getSubjectStart() {
	return (atoi(parts[subjectStartAt]));
}

unsigned CharSplitParser::getQueryEnd() {
	return (atoi(parts[queryEndAt]));
}

unsigned CharSplitParser::getSubjectEnd() {
	return (atoi(parts[subjectEndAt]));
}

unsigned CharSplitParser::getMatchArea() {
	return (atoi(parts[areaAt]));
}

float CharSplitParser::getPercentIdentity() {
	return (atof(parts[identityAt]));
}

bool CharSplitParser::getNextLine() {
	unsigned partCount;
	do {
		parts = puller->getCharSplitParts(splitter, partCount);
		if (parts == NULL) {
			// file data is exhausted so set partCount to exit loop
			partCount = minSplitCount;
		} else if (minSplitCount <= partCount) {
			// discard lines with identical query and subject
			if (0 == strcmp(parts[queryAt], parts[subjectAt])) {
				partCount = 0;
			}
		}
	} while (partCount < minSplitCount);
	return (parts != NULL);
}


CharSplitParser::~CharSplitParser() {
	delete(puller);
}
