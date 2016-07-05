/*
 * LineParser.cpp
 *
 *  Created on: Apr 28, 2010
 *      Author: mscott
 */

#include "LineParser.h"
#include <math.h>
#include <iostream>
#include <string.h>


LineParser::LineParser() {
}

float LineParser::getValueLog(float eValue, char *eText) {
	// check for exponential notation
	char *eSplit = strchr(eText, 'e');
	if (eSplit == NULL) {
		eSplit = strchr(eText, 'E');
	}
	if (eSplit == NULL) {
		double aid = atof(eText);
		// can not take log of 0.0 so current eValue will be return
		if (0.0 != aid) {
			aid = log10(aid);
			eValue = (float)aid;
		}
	} else {
		*eSplit = '\0';
		double aid = atof(eText);
		// can not take log of 0.0 so current eValue will be return
		if (aid != 0) {
			aid = log10(aid);
			++eSplit;
			// to adjust aid simply add exponent to return value
			aid += atof(eSplit);
			eValue = (float)aid;
		}
	}
	return eValue;
}

LineParser::~LineParser() {
}
