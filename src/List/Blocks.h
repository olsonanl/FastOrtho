/*
 * Blocks.h
 *
 *  Created on: Apr 8, 2010
 *      Author: mscott
 */

#ifndef BLOCKS_H_
#define BLOCKS_H_
#include "LinkFinder.h"
#include "BlockMarch.h"

class Blocks {
private:
	unsigned entrySize;
	unsigned blockSize;
	unsigned lastBlock;
	unsigned lastTop;
	void **blocksAt;
	unsigned atRoom;
	unsigned entriesInBlock;
	void* fullSort(void *catcher, void **inArray,
					int (*compare)(const char *left, const char *right));
	void* partSort(void *catcher, void **inArray, unsigned top,
					int (*compare)(const char *left, const char *right));
	void moveTail(BlockMarch *fromMarch, BlockMarch *storeMarch);
	void blockSort(void *catcher,
					int (*compare)(const char *left, const char *right));
	unsigned doubleRun(unsigned runBlocks,
			int (*compare)(const char *left, const char *right),
			BlockMarch *leftMarch, BlockMarch *rightMarch,
			BlockMarch *storeMarch, void **buildsAt);


public:
	Blocks(unsigned size);
	void add(char *);
	char* getEntry(long index);
	char* getAdd();
	void empty();
	void sort(int (*compare)(const char *left, const char *right));
	void keepPairs(
			bool (*merge)(const char *left, const char *right, char *saveIn));
	void demoteSingles(
			bool (*merge)(const char *left, const char *right, char *saveIn),
			Blocks *lesser, char *linkTrack);
	bool hasMembers();
	unsigned getEntryCount();
	BlockMarch* getFullMarch();
	char* getLastEntry();
	char* dropLastEntry(char *toDrop);
	char* emptyToOneBlock();
	bool chopBottom(LinkFinder* toFind);
	bool chopLower(LinkFinder* toFind, LinkFinder *belowPass, LinkFinder *above);
	bool chopTop(LinkFinder *below, LinkFinder* toFind);
	bool chopUpper(LinkFinder* toFind, LinkFinder *below, LinkFinder *abovePass);
	void absorb(Blocks *discard);



	virtual ~Blocks();
};

#endif /* BLOCKS_H_ */
