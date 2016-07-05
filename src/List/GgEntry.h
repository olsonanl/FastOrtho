/*
 * GgEntry.h
 *
 *  Created on: Apr 9, 2010
 *      Author: mscott
 */

#ifndef GGENTRY_H_
#define GGENTRY_H_

struct GgEntry {
public:
	char* name;
	unsigned parent;
	unsigned order;
};

#endif /* GGENTRY_H_ */
