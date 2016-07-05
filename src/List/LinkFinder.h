/*
 * LinkFinder.h
 *
 *  Created on: Apr 13, 2010
 *      Author: mscott
 */

#ifndef LINKFINDER_H_
#define LINKFINDER_H_

class LinkFinder {
public:
	unsigned low;
	unsigned high;

	bool found;
	bool aboveAll;
	bool belowAll;

	unsigned block;
	unsigned inBlock;
	unsigned entriesInBlock;
	char *blockStart;


	LinkFinder();

	void setToFind(unsigned left, unsigned right);

	void setLocation(unsigned blockSet, unsigned blockDepth,
				char* startSet, unsigned countSet);
	void setLocation(LinkFinder *model);

	char *getReader(unsigned entrySize);

	void copy(LinkFinder *toCopy);

	bool getMiddle(void **blocksAt, LinkFinder *above, LinkFinder *result);
	int compare(char *reader);


	virtual ~LinkFinder();
};

#endif /* LINKFINDER_H_ */
