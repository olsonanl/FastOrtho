/*
 * Links.h
 *
 *  Created on: Apr 9, 2010
 *      Author: mscott
 */

#ifndef LINKS_H_
#define LINKS_H_
#include "Blocks.h"
#include "../Parser/GgParser.h"
#include <vector>
#include "BlockMarch.h"


class Links {
private:
	unsigned *getParent;
	Blocks *storage;
	unsigned entrySize;
	unsigned *intraDepths;
	char **intraValues;

	unsigned taxonCount;

public:
	Links();
	Blocks *getStorage();
	void add(unsigned left, unsigned right, float rating);
	void add(char *toAdd);
	BlockMarch *getFullMarch();
	void absorb(Links *discard);
	void moveIntras(Links* bestInters, char *linkTrack);
	void sort();
	unsigned getEntryCount();
	void keepPairs();
	void demoteSingles(Links *lesser, char *linkTrack);
	void adjustForIntra(unsigned taxonCount, unsigned *toParent);
	void getClones(unsigned lowMatch, std::vector<unsigned> *result);
	virtual ~Links();
};

#endif /* LINKS_H_ */
