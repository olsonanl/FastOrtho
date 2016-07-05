/*
 * LinePull.h
 *
 *  Created on: Feb 26, 2010
 *      Author: mscott
 */

#ifndef LINEPULL_H_
#define LINEPULL_H_

#include <stdlib.h>
#include <stdio.h>
#include <string>

class LinePull {

private:
	char *buffer;
	char *pastBuffer;
	char *unused;
	size_t bufferRoom;
	FILE *reader;
	size_t available;

	bool didLastWhite;

	char **partPoints;
	unsigned partRoom;

	bool refreshBuffer();
	unsigned pullTries;

public:
	LinePull(size_t bufferRoom, const char *filePath);
	virtual ~LinePull();
	char *getLine();
	char ** getCharSplitParts(char splitter, unsigned &partCount);
	char *getFrontWhiteSplit();
	char *getNextWhiteSplit();
	char *trim(char *start, char *stop);
	char *trim(char *start);
};

#endif /* LINEPULL_H_ */
