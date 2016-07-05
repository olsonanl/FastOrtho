/*
 * LineParser.h
 *
 *  Created on: Apr 28, 2010
 *      Author: mscott
 */

#ifndef LINEPARSER_H_
#define LINEPARSER_H_
#include <vector>
#include "../FileIO/LinePull.h"
#define QUERY_AT 0
#define QUERY_START 1
#define QUERY_END 2
#define QUERY_LENGTH 3
#define SUBJECT_AT 4
#define SUBJECT_START 5
#define SUBJECT_END 6
#define SUBJECT_LENGTH 7
#define IDENTITY_AT 8
#define E_VALUE_AT 9
#define LENGTH_AT 10
#define MAP_COLUMN 11
#define INDEX_ROOM 12


class LineParser {
private:
	float minLogE;
public:
	LineParser();
	virtual char *getQueryName()= 0;
	virtual char *getSubjectName() = 0;
	virtual float getValueE() = 0;
	virtual char *getTextE() = 0;
	virtual unsigned getQueryStart() = 0;
	virtual unsigned getSubjectStart() = 0;
	virtual unsigned getQueryEnd() = 0;
	virtual unsigned getSubjectEnd() = 0;
	virtual unsigned getMatchArea() = 0;
	virtual float getPercentIdentity() = 0;
	virtual bool getNextLine() = 0;

	float getValueLog(float minLogE, char *eText);

	virtual ~LineParser();
};

#endif /* LINEPARSER_H_ */
