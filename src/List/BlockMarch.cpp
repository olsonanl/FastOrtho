/*
 * BlockMarch.cpp
 *
 *  Created on: Apr 15, 2010
 *      Author: mscott
 */

#include "BlockMarch.h"
#include <stdlib.h>
#include <iostream>
#include <string.h>

BlockMarch::BlockMarch(unsigned entrySize, unsigned blockSize,
					   unsigned atBlock, void **blockList,
					   unsigned lastAt, unsigned lastTop,
					   void **emptyCatch, unsigned *atAdd) {
	// save size of block entries
	stepSize = entrySize;
	// save size of block
	regularTop = blockSize;
	//  save index for current block
	at = atBlock;
	// save table of block locations
	jumpList = blockList;

	// lastAt is termination block
	//   and lastTop is termination position in lastAt
	if ((0 < lastAt) && (lastTop == 0)) {
		// desired termination is end of lastAt - 1
		--lastAt;
		lastTop = blockSize;
	}
	// mark termination block index
	finalAt = lastAt;
	//  mark termination position in last block
	finalTop = lastTop;

	// save means of moving emptied blocks to another Blocks object
	emptySave = emptyCatch;
	emptyAdd = atAdd;

	// assume there is nothing in march
	entry = NULL;
	if (at < finalAt) {
		// identify start of current block for possible emptySave transfer
		blockStart = jumpList[at];
		//  march begins at start of block
		entry = (char *)blockStart;
		// compute location when block change is required
		entryJump = entry + regularTop;
	} else if ((at == finalAt) && (0 < finalTop)) {
		// identify start of current block for possible emptySave transfer
		blockStart = jumpList[at];
		//  march begins at start of block
		entry = (char *)blockStart;
		// compute location that indicates end of march
		entryJump = entry + finalTop;
	}
}


void BlockMarch::reset(unsigned atBlock, unsigned lastAt, unsigned lastTop) {
	// set new current block
	at = atBlock;

	// if last block is empty try to retreat to full previous block
	if ((lastTop == 0) && (0 < lastAt)) {
		lastTop = regularTop;
		--lastAt;
	}

	//  set variables for detecting end of march
	finalAt = lastAt;
	finalTop = lastTop;
	// assume march is empty
	entry = NULL;
	if (at < finalAt) {
		// march starts in full block
		blockStart = jumpList[at];
		entry = (char *)blockStart;
		entryJump = entry + regularTop;
	} else if ((at == finalAt) && (0 < finalTop)) {
		// march starts in final block
		blockStart = jumpList[at];
		entry = (char *)blockStart;
		entryJump = entry + finalTop;
	}
}

unsigned BlockMarch::getAtBlockRoom() {
	// compute raw room in current block
	unsigned result = 0;
	if (entry != NULL) {
		result = entryJump - entry;
	}
	return result;
}

unsigned BlockMarch::getAtBlock() {
	return at;
}

unsigned BlockMarch::getAtFilled() {
	// get amount of block that has been passed over
	return (entry - (char *)blockStart);
}

char *BlockMarch::getEntry() {
	return entry;
}
char* BlockMarch::getNextEntry() {
	if (entry != NULL) {
		// advance to next entry in block
		entry += stepSize;
		if (entryJump <= entry) {
			// block is exhausted so switch block
			entry = advanceBlock();
		}
	}
	return entry;
}

void BlockMarch::transferTailBlocks(BlockMarch *sendTo) {
	// transfer is only supported if proper empty block transfer was
	//   specified in constructor
	if ((emptySave != NULL) && (emptySave == (sendTo->jumpList))) {
		// this should be called only after the remaining contents
		//    of startBlock have been used to fill the current
		//    block of sendTo

		// current block is considered consumed so code should
		//   transfer it to end of emptySave
		//   and position itself at next block
		advanceBlock();

		// skip over sendTo block that was filled earlier
		++(sendTo->at);

		// now complete block transfers are possible
		// compute number blocks to transfer
		unsigned tailBlocks = 1 + finalAt - at;
		//  transfer buffer blocks to location past transfer area
		if (0 < tailBlocks) {
			// full block address will be send to the jumpList of sendTo
			// so the empty blocks that are in the way of this transfer
			// need to be moved to the end of sendTo->jumpList
			unsigned emptyCount = *emptyAdd - (sendTo->at);

			// shift blocks at sendTo->at up tailBlocks in position
			memmove(emptySave + tailBlocks + (sendTo->at),
					emptySave + (sendTo->at), emptyCount * sizeof(void *));
			//  adjust indication of total blocks by amount that is to be
			//  transfered
			*emptyAdd += tailBlocks;
			// transfer blocks
			memcpy(emptySave + (sendTo->at), jumpList + at,
				   tailBlocks * sizeof(void *));
			at += tailBlocks;
			(sendTo->at) += tailBlocks;
		}
		sendTo->setBlockStart();
		// sendTo settings will only be accurate if this march
		//    completed in a full block but the all uses of this function
		//    are finished with sendTo in non-full block cases
	}
}


char* BlockMarch::advanceBlock() {
	// see if block should be transfered
	if (emptySave != NULL) {
		// move block
		emptySave[*emptyAdd] = blockStart;
		// adjust storage location for next emptied block
		++(*emptyAdd);
	}
	// advance to next block
	++at;
	// set per block variables
	setBlockStart();
	return entry;
}

void BlockMarch::setBlockStart() {
	if (finalAt < at) {
		// have consumed all blocks
		entry = NULL;
	} else {
		// get next block
		blockStart = jumpList[at];
		// set entry at start of next block
		entry = (char *)blockStart;
		// final block get different
		entryJump = entry;
		if (at < finalAt) {
			entryJump += regularTop;
		} else {
			entryJump += finalTop;
		}
	}
}


BlockMarch::~BlockMarch() {
	// TODO Auto-generated destructor stub
}
