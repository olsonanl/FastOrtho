/*
 * PickPiE.h
 *
 *  Created on: May 21, 2010
 *      Author: mscott
 */

#ifndef PICKPIE_H_
#define PICKPIE_H_

#include "LinkParser.h"
#include "../Parser/LineParser.h"

class PickPiE : public LinkParser {
private:
	LineParser *liner;

	bool initParts();
	bool setParts();
	void findNextIsland();
	void startNewQuery();
	bool findNextQuery();
	bool haveSameLeft();

	const char * lineLeft;
	const char * lineRight;

	unsigned matchArea;
	double piCut;

	bool sameLeft;
	bool haveParts;

public:
	PickPiE(LineParser *linerSet, GgParser *setGg,
			float setMaxE, double minPercent);
	bool setNextLine();
    bool isValid();
    bool hasBadPi();
	virtual ~PickPiE();
};

#endif /* PICKPIE_H_ */
