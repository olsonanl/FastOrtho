/*
 * LinkParser.h
 *
 *  Created on: Apr 28, 2010
 *      Author: mscott
 */

#ifndef LINKPARSER_H_
#define LINKPARSER_H_
#include "../Parser/GgParser.h"

class LinkParser {
public:
	LinkParser();
	virtual ~LinkParser();
	virtual bool setNextLine() = 0;
	virtual bool isValid() = 0;
	virtual bool hasBadPi() = 0;
	float minLogE;
	float maxLogE;
	bool atNewQuery;
	float eValue;
	double percentIdentity;
	double percentMatch;

	GgParser * ggParser;
	GgEntry *leftGene;
	GgEntry *rightGene;

};

#endif /* LINKPARSER_H_ */
