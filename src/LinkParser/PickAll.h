/*
 * PickAll.h
 *
 *  Created on: May 21, 2010
 *      Author: mscott
 */

#ifndef PICKALL_H_
#define PICKALL_H_

#include "LinkParser.h"
#include "../Parser/LineParser.h"

class PickAll : public LinkParser {
private:
	LineParser *liner;

	bool initParts();
	bool setParts();
	void findNextIsland();
	void startNewQuery();
	bool findNextQuery();
	bool haveSameLeft();
	void setMapSpan(unsigned start, unsigned end, unsigned *store);

	const char * lineLeft;
	const char * lineRight;

	Blocks *maps;

	bool sameLeft;
	bool haveParts;

	unsigned matchArea;
	double piCut;

	unsigned *geneLengths;

	unsigned queryLength;
	unsigned subjectLength;
	unsigned smallerLength;
	double matchCut;

public:
	PickAll(LineParser *linerSet, GgParser *setGg,
			float setMaxE, double minPercentIdentity, double minMatchCover);
	bool setNextLine();
	bool isValid();
	bool hasBadPi();
	virtual ~PickAll();
};

#endif /* PICKALL_H_ */
