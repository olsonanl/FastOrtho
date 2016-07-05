/*
 * BpoAllParser.h
 *
 *  Created on: May 21, 2010
 *      Author: mscott
 */

#ifndef BPOALLPARSER_H_
#define BPOALLPARSER_H_

#include "GgParser.h"
#include "LineParser.h"
#include "../FileIO/LinePull.h"

class BpoAllParser : public LineParser {
private:
	LinePull *puller;
	char splitter;
	char **parts;
	char *nextMap;
	unsigned queryAt;
	unsigned queryStart;
	unsigned queryEnd;
	unsigned queryLength;
	unsigned subjectAt;
	unsigned subjectStart;
	unsigned subjectEnd;
	unsigned subjectLength;
	unsigned eValueAt;
	unsigned identityAt;
	unsigned mapColumn;
	unsigned minColumnCount;
	bool getMapValues(char *mapAt);
	void setIndex(unsigned findAt, unsigned *settings,
							unsigned *toSet);
	float eReturn;
	float minLogE;

	GgParser *ggParser;
	const char *lastQueryName;
	const char *lastSubjectName;

public:
	BpoAllParser(GgParser* ggSet, const char *bpoFile, char splitSet,
				 double maxWeight, unsigned *indexSettings);
	char *getQueryName();
	char *getSubjectName();
	float getValueE();
	char* getTextE();
	unsigned getMatchArea();
	unsigned getQueryStart();
	unsigned getSubjectStart();
	unsigned getQueryEnd();
	unsigned getSubjectEnd();
	float getPercentIdentity();
	bool getNextLine();

	virtual ~BpoAllParser();

};

#endif /* BPOALLPARSER_H_ */
