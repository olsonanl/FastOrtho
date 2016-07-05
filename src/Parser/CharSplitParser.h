/*
 * CharSplitParser.h
 *
 *  Created on: Apr 28, 2010
 *      Author: mscott
 */

#ifndef CHARSPLITPARSER_H_
#define CHARSPLITPARSER_H_

#include "LineParser.h"
#include "../FileIO/LinePull.h"

class CharSplitParser : public LineParser {
private:
	LinePull *puller;
	char splitter;
	char **parts;
	float minLogE;
	unsigned minSplitCount;
	unsigned queryAt;
	unsigned subjectAt;
	unsigned eValueAt;
	unsigned identityAt;
	unsigned areaAt;
	unsigned queryStartAt;
	unsigned queryEndAt;
	unsigned subjectStartAt;
	unsigned subjectEndAt;
	void setIndex(unsigned findAt, unsigned *settings, unsigned *toSet);

public:
	CharSplitParser(char splitSet, const char *splitFile, float minLog,
					unsigned *indexSettings);

	char *getQueryName();
	char *getSubjectName();
	float getValueE();
	char *getTextE();
	float getPercentIdentity();
	unsigned getMatchArea();
	unsigned getQueryStart();
	unsigned getSubjectStart();
	unsigned getQueryEnd();
	unsigned getSubjectEnd();
	bool getNextLine();


	virtual ~CharSplitParser();
};

#endif /* CHARSPLITPARSER_H_ */
