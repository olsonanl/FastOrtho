/*
 * PickE.h
 *
 *  Created on: May 19, 2010
 *      Author: mscott
 */

#ifndef PICKE_H_
#define PICKE_H_
#include "LinkParser.h"
#include "../Parser/LineParser.h"

class PickE : public LinkParser {
private:
	LineParser *liner;

	bool initParts();
	bool setParts();
	void startNewQuery();
	void checkExtract();
	void findNextIsland();
	bool haveSameLeft();

	const char * lineLeft;
	const char * lineRight;

	bool sameLeft;
	bool haveParts;

public:
	PickE(LineParser *linerSet, GgParser *setGg, float setMaxE);
	bool setNextLine();
	bool isValid();
	bool hasBadPi();
	bool findNextQuery();

	virtual ~PickE();
};

#endif /* PICKE_H_ */
