/*
 * BlockMarch.h
 *
 *  Created on: Apr 15, 2010
 *      Author: mscott
 */

#ifndef BLOCKMARCH_H_
#define BLOCKMARCH_H_

class BlockMarch {
private:
	// address of current block
	void *blockStart;
	//  entry position in current block
	char *entry;
	//  termination location of current block
	char *entryJump;
	//  room of entry in current block (entry += stepSize to advance)
	unsigned stepSize;

	// list of blocks locations
	void **jumpList;

	// index off current block blockStart = jumpList[at]
	unsigned at;

	// index of last usable block
	unsigned finalAt;
	// termination location in jumpList[finalAt]
	unsigned atStop;

	// usable room jumpList[finalAt]
	unsigned finalTop;
	// usable room for other blocks in jumpList
	unsigned regularTop;

	// equivalent of blockList but in another Blocks object
	void **emptySave;
	//  provides position in emptySave for transfer of blockList members
	//    that are emptied by march
	unsigned *emptyAdd;

public:
	BlockMarch(unsigned entrySize, unsigned blockSize,
			   unsigned atBlock, void **blockList,
			   unsigned lastAt, unsigned lastTop,
			   void **emptyCatch, unsigned *atAdd);
	char *getEntry();
	char *getNextEntry();
	char *advanceBlock();
	void setBlockStart();
	void transferTailBlocks(BlockMarch *);
	unsigned getAtBlock();
	unsigned getAtFilled();

	void reset(unsigned atBlock, unsigned lastAt, unsigned lastTop);
	unsigned getAtBlockRoom();
	virtual ~BlockMarch();
};

#endif /* BLOCKMARCH_H_ */
