/*
 * LinkFinder.cpp
 *
 *  Created on: Apr 13, 2010
 *      Author: mscott
 */

#include "LinkFinder.h"
#include <stdlib.h>
#include <stdio.h>
#include <iostream>

LinkFinder::LinkFinder() {
}

void LinkFinder::setToFind(unsigned left, unsigned right) {
	// set coordinates of desired link
	low = left;
	high = right;
}

void LinkFinder::setLocation(unsigned blockSet, unsigned blockDepth,
							 char* startSet, unsigned countSet) {
	// search staring point for search relative to some Blocks object
	block = blockSet;  // index of block
	inBlock = blockDepth; // entry depth in block
	entriesInBlock = countSet;  // number of entries in block
	blockStart = startSet;  // location of block
}

void LinkFinder::setLocation(LinkFinder *model) {
	block = model->block;
	inBlock = model->inBlock;
	entriesInBlock = model->entriesInBlock;
	blockStart = model->blockStart;
}

char* LinkFinder::getReader(unsigned entrySize) {
	char *reader = blockStart;
	unsigned add = entrySize * inBlock;
	reader += add;
	return (blockStart + (entrySize * inBlock));
}

bool LinkFinder::getMiddle(void **blocksAt, LinkFinder *above,
							LinkFinder *result) {
	// assume that no entries are located between this and above
	bool isSet = false;
	if (block == above->block) {
		// both this & above point entries in the same block
		//  get the mean of their entry indices in their common block
		unsigned depthSet = ((inBlock + above->inBlock) >> 1);
		// if mean entry index matches this index then
		//   there is no middle
		if (depthSet != inBlock) {
			// set return with middle location
			result->setLocation(block, depthSet, blockStart, entriesInBlock);
			//  indicate that a middle was found
			isSet = true;
		}
	} else {
		long midFind = (above->block) - block;
		midFind *= entriesInBlock;
		// midFind = # of entries from start of this block
		//    to start of above block

		// drop front of this block
		midFind -= inBlock;
		//  add front of above block
		midFind += above->inBlock;
		midFind >>= 1;
		// midFind is the number of entries to add to this
		if (0 < midFind) {
			midFind += inBlock;
			// midFind is the number of entries past the start of
			//    the front of this block for the middle
			if (midFind < (long)entriesInBlock) {
				// middle is in the same block as this
				result->setLocation(block, (unsigned)midFind, blockStart,
									entriesInBlock);
			} else {
				// divide to get number of blocks past ths
				//   where we will find the entry associated with midFind
				//   remainder of division determines entry midFind
				//   index in its block
				unsigned blockSet = (unsigned)(midFind / entriesInBlock);
				unsigned depthSet = (unsigned)(midFind % entriesInBlock);
				blockSet += block;
				// split to allow potential different entry size
				//   if midFind ends up in final block of Blocks object
				if (blockSet < above->block) {
					result->setLocation(blockSet, depthSet,
								(char *)(blocksAt[blockSet]),
								entriesInBlock);
				} else {
					result->setLocation(blockSet, depthSet, above->blockStart,
										above->entriesInBlock);
				}
			}
			isSet = true;
		}
	}
	return isSet;
}



int LinkFinder::compare(char *values) {
	unsigned check = ((unsigned *)values)[0];
	int result = 0;
	if (low < check) {
		--result;
	} else if (check < low) {
		++result;
	} else {
		check = ((unsigned *)values)[1];
		if (high < check) {
			--result;
		} else if (check < high) {
			++result;
		}
	}
	return result;
}

LinkFinder::~LinkFinder() {
}
