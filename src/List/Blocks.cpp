/*
 * Blocks.cpp
 *
 *  Created on: Apr 8, 2010
 *      Author: mscott
 */

#include "Blocks.h"
#include "BlockMarch.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <iostream>

Blocks::Blocks(unsigned size) {
	// set last indicators for a Blocks without entries
	lastBlock = 0;
	lastTop = 0;
	// compute block size based on size
	entrySize = size;
	entriesInBlock = 0x10000;
	blockSize = size * entriesInBlock;
	//  create array for tracking blocks used
	atRoom = 0x100;
	blocksAt = (void **)malloc(atRoom * sizeof(void *));
	// allocate memory for first block
	*blocksAt = malloc(blockSize);
}

bool Blocks::hasMembers() {
	return (0 < lastTop);
}

unsigned Blocks::getEntryCount() {
	unsigned result = lastBlock * entriesInBlock;
	result += (lastTop / entrySize);
	return result;
}

void Blocks::add(char *entry) {
	if (lastTop < blockSize) {
		char* storeAt = (char *)(*(blocksAt + lastBlock));
		storeAt += lastTop;
		memcpy(storeAt, entry, entrySize);
		lastTop += entrySize;
	} else {
		++lastBlock;
		if (atRoom <= lastBlock) {
			atRoom += 0x100;
			blocksAt = (void **)realloc(blocksAt, atRoom * sizeof(void *));
		}
		void * newBlock = malloc(blockSize);
		*(blocksAt + lastBlock) = newBlock;
		memcpy((char *)newBlock, entry, entrySize);
		lastTop = entrySize;
	}
}

char* Blocks::emptyToOneBlock() {
	char *result = NULL;
	// keep NULL return if there are no entries
	if (0 < lastTop) {
		unsigned entryCount = getEntryCount();
		unsigned room = entryCount * entrySize;

		// create a single section of ram that is large enough
		//   to hold all entries
		result = (char *)(malloc(room));
		if (result == NULL) {
			printf("one block to %d malloc failure", room);
			std::cout << std::endl;
		}

		// move all of the full blocks to result
		char *storeAt = result;
		char *from;
		for (unsigned i = 0; i < lastBlock; i++) {
			from = (char *)(blocksAt[i]);
			memcpy(storeAt, from, blockSize);
			storeAt += blockSize;
			// release ram used for full block
			free(from);
		}
		// move last block to result
		from = (char *)(blocksAt[lastBlock]);
		memcpy(storeAt, from, lastTop);

		// alway keep on allocated block in object
		*blocksAt = (void *)from;
		// indicate that object has no entries
		lastBlock = 0;
		lastTop = 0;
	}
	return result;
}


char* Blocks::getLastEntry() {
	char *result = NULL;
	if (0 < lastTop) {
		result = (char *)(*(blocksAt + lastBlock));
		result += lastTop - entrySize;
	}
	return result;
}

char* Blocks::dropLastEntry(char *toDrop) {
	char *result = NULL;
	if (0 < lastTop) {
		lastTop -= entrySize;
		if (0 < lastTop) {
			// result is just before toDrop in same block
			result = toDrop - entrySize;
		} else {
			// toDrop was the last member in its block
			if (0 < lastBlock) {
				// release block which held toDrop
				free(toDrop);
				// move lastBlock to earlier member
				--lastBlock;
				lastTop = blockSize;
				// point result at last entry in new lastBlock
				result = (char *)(blocksAt[lastBlock]);
				result += blockSize - entrySize;
			}
		}
	}
	return result;
}

void Blocks::empty() {
	for (unsigned i = 1; i <= lastBlock; i++) {
		free(blocksAt + i);
	}
	lastBlock = 0;
	lastTop = 0;
}

char* Blocks::getAdd() {
	char *result;
	if (lastTop < blockSize) {
		// can store in current block
		result = (char *)(blocksAt[lastBlock]);
		result += lastTop;
		// set for entry caller is about to set
		lastTop += entrySize;
	} else {
		// next store will go in a new block
		++lastBlock;
		if (atRoom <= lastBlock) {
			// need more room for tracking all of the blocks
			atRoom += 0x100;
			blocksAt = (void **)realloc(blocksAt, atRoom * sizeof(void *));
		}
		void *newBlock = malloc(blockSize);
		blocksAt[lastBlock] = newBlock;
		result = (char *)newBlock;
		// set for entry caller is about to set
		lastTop = entrySize;
	}
	return result;
}


void Blocks::keepPairs(
		bool(*merge)(const char *left, const char *right, char *saveIn)) {
	if (lastBlock == 0) {
		// no need to worry about block changing
		char *preAt = (char *)(*blocksAt);
		char *saveIn = preAt;
		char *postAt = preAt + entrySize;
		// mark end of block
		char *end = preAt + lastTop;
		while (postAt < end) {
			if (merge(preAt, postAt, saveIn)) {
				// left and right have been merged
				//   to a single entry at saveIn

				// prepare saveIn for next merge
				saveIn += entrySize;
				//  advance postAt so next
				//  pair comparison will be
				//    after postAt
				postAt += entrySize;
			}
			//  shift pair to compare
			preAt = postAt;
			postAt += entrySize;
		}
		//  adjust so only merged pairs are kept
		lastTop = saveIn - (char *)(*blocksAt);
	} else {
		char *preAt = (char *)(*blocksAt);
		char *saveIn = preAt;
		char *saveEnd = preAt + blockSize;
		char *end = saveEnd;
		unsigned saveBlock = 0;
		char *postAt = preAt + entrySize;
		unsigned postBlock = 0;

		// lack of pairs to compare is
		//     signaled by lastBlock < postBlock
		while (postBlock <= lastBlock) {
			if ((merge(preAt, postAt, saveIn))) {
				//  saveIn contains merged version
				//   of preAt & postAt

				//  advance saveIn to collect next
				//   merge pair
				saveIn += entrySize;
				if (saveEnd <= saveIn) {
					// advance block
					++saveBlock;
					// set saveIn at end of next block
					saveIn = (char *)(*(blocksAt + saveBlock));
					//  will never reach end of last block
					//    so only need to set saveEnd for full blocks
					saveEnd = saveIn;
					saveEnd += blockSize;
				}
				// advance postAt
				postAt += entrySize;
				if (end <= postAt) {
					++postBlock;
					// check to see if postAt is at last entry
					if (postBlock <= lastBlock) {
						// advance postAt to next block
						postAt = (char *)(blocksAt [postBlock]);
						//  set end of postBlock
						end = postAt;
						if (postBlock < lastBlock) {
							end += blockSize;
						} else {
							end += lastTop;
						}
					}
				}
			}
			// without merge new pre is old post
			// in merge code post was advanced
			//      so post value still represents next pre
			preAt = postAt;

			// advance post to pair with pre
			postAt += entrySize;
			if (end <= postAt) {
				// need to advance post to start of next block
				++postBlock;
				if (postBlock <= lastBlock) {
					// set postAt at start of next block
					postAt = (char *)(blocksAt [postBlock]);
					// mark end of block
					end = postAt;
					if (postBlock < lastBlock) {
						end += blockSize;
					} else {
						end += lastTop;
					}
				}
			}
		}

		// have finished merge
		if (saveBlock == 0) {
			// set lastTop
			lastTop = saveIn - (char *)(*blocksAt);
		} else {
			// compute the amount in merged into saveBlock
			char *blockStart = (char *)(blocksAt[saveBlock]);
			lastTop = saveIn - blockStart;
			// if block is empty retreat to previous block
			if (lastTop == 0) {
				lastTop = blockSize;
				--saveBlock;
			}
		}

		//  release blocks which did not get any merge entries
		unsigned release = saveBlock + 1;
		while (release <= lastBlock) {
			free(blocksAt[release]);
			++release;
		}
		// adjust lastBlock to agree with end of merged contents
		lastBlock = saveBlock;
	}
}

void Blocks::demoteSingles(
		bool(*merge)(const char *left, const char *right, char *saveIn),
		Blocks *lesser, char *linkTrack) {
	BlockMarch* postMarch = getFullMarch();
	BlockMarch* saveMarch = getFullMarch();
	char* preAt = postMarch->getEntry();
	char* saveAt = preAt;
	char *postAt = postMarch->getNextEntry();
	while (postAt != NULL) {
		// compare consecutive entries
		//    and merge into a single entry if
		//      they link the same 2 genomes
		if (merge(preAt, postAt, saveAt)) {
			// record that both genomes have best inter links
			unsigned index = ((unsigned *)saveAt)[0];
			linkTrack[index] = '\1';
			index = ((unsigned *)saveAt)[1];
			linkTrack[index] = '\1';
			// advance saveAt for collecting next merge
			saveAt = saveMarch->getNextEntry();
			//  (preAt, postAt) got merged so
			//    next pair will start at postAt
			preAt = postMarch->getNextEntry();
			if (preAt == NULL) {
				// all entries are gone so set for exit
				postAt = NULL;
			} else {
				// set postAt just past preAt
				postAt = postMarch->getNextEntry();
			}
		} else {
			//  preAt is a singleton so save in lesser
			lesser->add(preAt);
			// next pair starts at non-matching postAt
			preAt = postAt;
			// get next entry to complete pair
			postAt = postMarch->getNextEntry();
		}
	}
	// check to see if last entry was a singleton
	if (preAt != NULL) {
		lesser->add(preAt);
	}
	// release postMarch resources
	delete (postMarch);

	// adjust this object to only hold pair merges
	lastTop = saveMarch->getAtFilled();
	unsigned saveBlock = saveMarch->getAtBlock();
	// release blocks that did not get any merged entries
	for (unsigned i = saveBlock + 1; i <= lastBlock; i++) {
		free(blocksAt[i]);
	}
	// set new last block after release
	lastBlock = saveBlock;
	// release saveMarch resources
	delete (saveMarch);
}

BlockMarch* Blocks::getFullMarch() {
	return (new BlockMarch(entrySize, blockSize, 0, blocksAt,
							lastBlock, lastTop, NULL, NULL));
}

void Blocks::sort(int (*compare)(const char *left, const char *right)) {
	void *catcher = malloc(blockSize);
	if (lastBlock == 0) {
		// only 1 block
		catcher = partSort(catcher, blocksAt, lastTop, compare);
		// release extra block used for merging
		free(catcher);
	} else {
		// make sure all blocks are sorted
		for (unsigned i = 0; i < lastBlock; i++) {
			catcher = fullSort(catcher, blocksAt + i, compare);
		}
		catcher = partSort(catcher, blocksAt + lastBlock, lastTop, compare);
		// merge whole blocks
		blockSort(catcher, compare);
	}
}

void Blocks::blockSort(void *catcher,
						int (*compare)(const char *left, const char *right))
{
	// merging may require up to 2 more blocks than
	//  are currently being used
	void **buildsAt = (void **)malloc((3 + lastBlock) * sizeof(void *));

	//  set start of merge with unused input block
	*buildsAt = catcher;
	//  add an additional unused block to merge area
	*(buildsAt + 1) = malloc(blockSize);

	// set to add more blocks to buildsAt as they are emptied in blocksAt
	unsigned buildAdd = 2;

	//  variable for number of blocks to merge
	// convenience for checking block overflows
	// set controls for storing into buildsAt
	unsigned secondTop = lastTop;
	if (1 < lastBlock) {
		secondTop = blockSize;
	}
	// set storeMarch to collect in merge area
	BlockMarch storeMarch(entrySize, blockSize, 0, buildsAt,
						 lastBlock, lastTop, NULL, NULL);
	// set leftMarch to pull from 0 block
	BlockMarch leftMarch(entrySize, blockSize, 0, blocksAt, 0, blockSize,
						 buildsAt, &buildAdd);

	// set rightMarch to pull from 1 block
	BlockMarch rightMarch(entrySize, blockSize, 1, blocksAt, 1, secondTop,
						  buildsAt, &buildAdd);
	// merge first to blocks
	unsigned runBlocks =
		doubleRun(1, compare, &leftMarch, &rightMarch, &storeMarch,
					buildsAt);
	while (runBlocks <= lastBlock) {
		// doubleRun will finish with 2 empty blocks at start of buildsAt
		buildAdd = 2;
		//  set leftMarch to span 1st sorted section
		leftMarch.reset(0, runBlocks - 1, blockSize);
		// find block just past second sorted section
		unsigned rightStop = (runBlocks << 1);
		if (rightStop <= lastBlock) {
			// right section is composed of full blocks
			rightMarch.reset(runBlocks, rightStop - 1, blockSize);
		} else {
			//  right section goes to end of area to sort
			rightMarch.reset(runBlocks, lastBlock, lastTop);
		}
		// storeMarch will be full when its size matches current
		//   size of this object
		storeMarch.reset(0, lastBlock, lastTop);

		runBlocks =
			doubleRun(runBlocks, compare, &leftMarch, &rightMarch,
							&storeMarch, buildsAt);
	}
	// release memory in extra blocks
	free(buildsAt[0]);
	free(buildsAt[1]);
	// release merge block tracking
	free(buildsAt);
}

unsigned Blocks::doubleRun(unsigned runBlocks,
		int (*compare)(const char *left, const char *right),
		BlockMarch *leftMarch, BlockMarch *rightMarch,
		BlockMarch *storeMarch, void **buildsAt) {
	unsigned leftBlock = 0;
	while (leftBlock <= lastBlock) {
		char *storeAt = storeMarch->getEntry();
		char* leftAt = leftMarch->getEntry();
		char* rightAt = rightMarch->getEntry();
		// loop until merging is finished
		while ((leftAt != NULL) && (rightAt != NULL)) {
			// check to see which gets merged next
			int differ = compare(leftAt, rightAt);
			if (differ <= 0) {
				// merge left entry
				memcpy(storeAt, leftAt, entrySize);
				leftAt = leftMarch->getNextEntry();
				storeAt = storeMarch->getNextEntry();
			}
			if (0 <= differ) {
				memcpy(storeAt, rightAt, entrySize);
				rightAt = rightMarch->getNextEntry();
				storeAt = storeMarch->getNextEntry();
			}
		}
		if (leftAt != NULL) {
			moveTail(leftMarch, storeMarch);
		} else if (rightAt != NULL) {
			moveTail(rightMarch, storeMarch);
		}
		// skip over merged sections
		leftBlock += (runBlocks << 1);
		if (leftBlock <= lastBlock) {
			// have more to merge
			unsigned rightBlock = leftBlock + runBlocks;
			if (lastBlock < rightBlock) {
				//  just tack on sorted tail
				leftMarch->reset(leftBlock, lastBlock, lastTop);
				// move entire leftMarch area to merge
				moveTail(leftMarch, storeMarch);
				// set leftBlock to exit loop
				leftBlock = lastBlock + 1;
			} else {
				// set next left section for merging
				leftMarch->reset(leftBlock, rightBlock - 1, blockSize);
				//  compute end of full right merge
				unsigned rightStop = rightBlock + runBlocks - 1;

				if (rightStop < lastBlock) {
					// have enough left for a full right block
					rightMarch->reset(rightBlock, rightStop, blockSize);
				} else {
					//  set rightMarch to merge tail of from area
					rightMarch->reset(rightBlock, lastBlock, lastTop);
				}
			}
		}
	}
	// move merge blocks into blocksAt for this object
	memcpy(blocksAt, buildsAt, (lastBlock + 1) * (sizeof(void *)));
	// move extra blocks to lowest positions in buildsAt
	memcpy(buildsAt, (buildsAt + lastBlock + 1), 2 * sizeof(void *) );

	runBlocks <<= 1;
	return runBlocks;
}

void Blocks::moveTail(BlockMarch *fromMarch, BlockMarch *storeMarch) {
	// compute amount remaining in from and store blocks
	unsigned fromSize = fromMarch->getAtBlockRoom();
	unsigned storeSize = storeMarch->getAtBlockRoom();
	char *storeAt = storeMarch->getEntry();
	char *fromAt = fromMarch->getEntry();
	do {
		if (fromSize < storeSize) {
			// store has enough room to empty from block
			memcpy(storeAt, fromAt, fromSize);
			// adjust storeAt to reflect copy
			storeAt += fromSize;
			//  compute amount left in storeBlock
			storeSize -= fromSize;
			//  set fromAt at start of next from block
			//     and adjust its current block
			fromAt = fromMarch->advanceBlock();
			fromSize = fromMarch->getAtBlockRoom();

		} else if (fromSize == storeSize) {
			// from block will exactly fill store block
			//    any remaining blocks can be handled
			//    by adjusting buildsAt instead of moving
			//    actual entries
			memcpy(storeAt, fromAt, fromSize);
			fromMarch->transferTailBlocks(storeMarch);
			// set to exit loop
			fromSize = 0;
		} else {
			// fill remainder of store block
			memcpy(storeAt, fromAt, storeSize);
			// adjust fromAt in response tof move
			fromAt += storeSize;
			fromSize -= storeSize;
			//  move storeAt to next block and adjust its current block
			storeAt = storeMarch->advanceBlock();
			// set storeSize with size of next block
			storeSize = storeMarch->getAtBlockRoom();
		}
	} while (0 < fromSize);
}



void* Blocks::fullSort(void *catcher, void **inArray,
		int (*compare)(const char *left, const char *right)) {
	// entries have yet to be sorted
	unsigned runSize = entrySize;
	while (runSize < blockSize) {
		// access data with runSize sorted entry sections
		void *toSort = *inArray;
		//  set to other block for catching merge
		char *storeAt = (char *)catcher;
		//  set stop to mark filling of merge block
		char *storeStop = storeAt + blockSize;
		//  left at left of source block
		char *leftAt = (char *)toSort;
		while (storeAt < storeStop) {
			//  mark end of stored section
			char *leftStop = leftAt + runSize;
			//  mark ends of other sorted section
			char *rightAt = leftStop;
			char *rightStop = rightAt + runSize;
			int differ;
			// loop until at least one sorted section is merged
			while ((leftAt < leftStop) && (rightAt < rightStop)) {
				// determine which should be placed into merge next
				differ = compare(leftAt, rightAt);
				if (differ <= 0) {
					memcpy(storeAt, leftAt, entrySize);
					leftAt += entrySize;
					storeAt += entrySize;
				}
				if (0 <= differ) {
					memcpy(storeAt, rightAt, entrySize);
					rightAt += entrySize;
					storeAt += entrySize;
				}
			}
			if (leftAt < leftStop) {
				// send remainder of left section into merge
				differ = leftStop - leftAt;
				memcpy(storeAt, leftAt, differ);
				storeAt += differ;
			} else if (rightAt < rightStop) {
				// send remainder of right section into merge
				differ = rightStop - rightAt;
				memcpy(storeAt, rightAt, differ);
				storeAt += differ;
			}
			leftAt = rightStop;
		}
		// replace entry in blocksAt with
		//    block with larger runSize
		*inArray = catcher;
		//  set less sorted block to catch next merge
		catcher = toSort;
		runSize <<= 1;
	}
	return catcher;
}



void* Blocks::partSort(void *catcher, void **inArray, unsigned top,
		int (*compare)(const char *left, const char *right)) {
	// sorted runs start at 1 entry
	unsigned runSize = entrySize;
	while (runSize < top) {
		// get area to sort
		void *toSort = *inArray;
		// place version with double sorted lengths in catcher
		char *storeAt = (char *)catcher;
		//  set marker for testing when catcher is full
		char *storeStop = storeAt + top;
		//  set left at leftmost position
		char *leftAt = (char *)toSort;
		// set marker for testing when toSort is empty
		char *end = leftAt + top;
		// loop until catcher is full
		while (storeAt < storeStop) {
			// put leftStop at end of sorted section
			char *leftStop = leftAt + runSize;
			if (end <= leftStop) {
				// only left section so complete fill of catcher
				memcpy(storeAt, leftAt, end - leftAt);
				// indicate that catcher is full
				storeAt = storeStop;
			} else {
				// set right section for merging
				char *rightAt = leftStop;
				char *rightStop = rightAt + runSize;
				// make sure right section does not overflow data
				if (end < rightStop) {
					rightStop = end;
				}
				int differ;
				// loop until at least one merging section is emptied
				while ((leftAt < leftStop) && (rightAt < rightStop)) {
					differ = compare(leftAt, rightAt);
					if (differ <= 0) {
						// place left entry in merge
						memcpy(storeAt, leftAt, entrySize);
						leftAt += entrySize;
						storeAt += entrySize;
					}
					if (0 <= differ) {
						// place right entry in merge
						memcpy(storeAt, rightAt, entrySize);
						rightAt += entrySize;
						storeAt += entrySize;
					}
				}
				if (leftAt < leftStop) {
					// finish end of merge with remainder of left
					differ = leftStop - leftAt;
					memcpy(storeAt, leftAt, differ);
					storeAt += differ;
				} else if (rightAt < rightStop) {
					// finish end of merge with remainder of right
					differ = rightStop - rightAt;
					memcpy(storeAt, rightAt, differ);
					storeAt += differ;
				}
				// set left just after emptied right section
				leftAt = rightStop;
			}
		}
		// update array with version that has longer sorted sections
		*inArray = catcher;
		//  set catcher with block that was just emptied
		catcher = toSort;
		// adjust runSize to show consequence of section merging
		runSize <<= 1;
	}
	return catcher;
}

char* Blocks::getEntry(long index) {
	unsigned depth = index & 0xFFFF;
	depth *= entrySize;
	index >>= 16;
	char *result = (char *)(*(blocksAt + index));
	result += depth;
	return result;
}

// called with smallest desired target
//   assumes that the Blocks entries are sorted in ascending order
//   if toFind is not found among the Blocks entries the
//      location of toFind will describe the smallest Blocks entry
//      that is larger than toFind
bool Blocks::chopBottom(LinkFinder* toFind) {
	toFind->found = false;
	if (lastTop == 0) {
		// no entries so nothing to find
		toFind->aboveAll = true;
		toFind->belowAll = true;
	} else {
		// get blockStart for 0 entry
		char *smallest = (char *)(blocksAt[0]);
		// compare smallest to toFind
		int differ = toFind->compare(smallest);
		if (differ < 0) {
			//  toFind is smaller than the entire block
			//   so the smallest block member is the desired return
			toFind->blockStart = smallest;
			toFind->block = 0;
			toFind->inBlock = 0;
			if (lastBlock == 0) {
				toFind->entriesInBlock = lastTop / entrySize;
			} else {
				toFind->entriesInBlock = entriesInBlock;
			}
			// to find is below entire block
			toFind->belowAll = true;
			toFind->aboveAll = false;
		} else if (differ == 0) {
			// no more work is required after a find
			toFind->found = true;
		} else {
			//  smallest < toFind
			toFind->belowAll = false;
			// create bounds to support binary search logic
			bool doLoop = false;
			LinkFinder below;
			LinkFinder above;
			if (lastBlock == 0) {
				// only 1 block to consider

				// get access to largest block member
				// single block in this means its block start = smallest
				char *reader = smallest;
				reader += lastTop - entrySize;
				differ = toFind->compare(reader);
				if (0 < differ) {
					// toFind is larger than entire block
					//    so use aboveAll to indicate failure
					//    and chop entire block
					toFind->aboveAll = true;
				} else if (differ == 0) {
					//  no more work if a match is found
					toFind->found = true;
				} else {
					//  set bounds for binary search
					unsigned entryCount = lastTop / entrySize;
					below.setLocation(0, 0, smallest, entryCount);
					above.setLocation(0, entryCount - 1, smallest, entryCount);
					doLoop = true;
				}
			} else {
				// get blockStart for last block
				char *biggest = (char *)(blocksAt[lastBlock]);
				char *reader = biggest;
				reader += lastTop - entrySize;
				// compare toFind to largest
				differ = toFind->compare(reader);
				if (0 < differ) {
					//  largest < toFind so nothing to return
					toFind->aboveAll = true;
				} else if (differ == 0) {
					// no more work if a match is found
					toFind->found = true;
				} else {
					// set bounds for middle search
					unsigned entryCount = lastTop / entrySize;
					below.setLocation(0, 0, smallest, entriesInBlock);
					above.setLocation(lastBlock, entryCount - 1, biggest,
									  entryCount);
					doLoop = true;
				}
			}
			//  see if binary search is required
			if (doLoop) {
				// since above is set toFind is not larger than all entries
				toFind->aboveAll = false;
				//  get middle for testing
				LinkFinder toCheck;
				doLoop = below.getMiddle(blocksAt, &above, &toCheck);
				while (doLoop) {
					char *reader = toCheck.getReader(entrySize);
					differ = toFind->compare(reader);
					if (differ < 0) {
						//  have a new lower bound for above
						above.setLocation(&toCheck);
						// continue binary search logic
						doLoop =
							below.getMiddle(blocksAt, &above, &toCheck);
					} else if (differ == 0) {
						//  match ends work
						toFind->found = true;
						doLoop = false;
					} else {
						// have better bound for below
						below.setLocation(&toCheck);
						// continue binary search logic
						doLoop =
							below.getMiddle(blocksAt, &above, &toCheck);
					}
				}
				if (!(toFind->found)) {
					toFind->setLocation(&above);
				}
			}
		}
	}
	return (toFind->found);
}

//  toFind describes the smallest desired target that has not been
//     eliminated by a search failure
//  the below to aboveModel range will be check to see if it contains
//     an entry that matches toFind.
//  If no matching entry exists the potential exist to reduce this
//     range for future searches by adjusting the location of below
//  true is returned if a match is found for toFind
bool Blocks::chopLower(LinkFinder* toFind, LinkFinder* below,
					   LinkFinder* aboveModel) {

	//  aboveModel should not be modified
	//    so create a copy to support binary reduction of search bounds
	LinkFinder above;
	above.setLocation(aboveModel);

	// get access to initial below
	char *reader = below->getReader(entrySize);
	// compare smallest bound to toFind
	int differ = toFind->compare(reader);
	toFind->aboveAll = false;
	//  below might also match toFind and set result
	toFind->found = (differ == 0);
	// differ == 0 means search is complete and no more work is required
	//  differ < 0 means toFind is smaller than all entries in
	//    our range so below should be returned unmodified
	if (0 < differ) {
		// toFind will reduce bounds
		char *reader = above.getReader(entrySize);
		differ = toFind->compare(reader);
		if (0 < differ) {
			// entire block should be discarded
			//   an no more targets need be checked
			toFind->aboveAll = true;
		} else if (differ == 0) {
			toFind->found = true;
		} else {
			//  the below entry is smaller than toFind
			//     and the above entry is larger than toFind
			//   it is time to try a binary search
			LinkFinder toCheck;
			bool doLoop = below->getMiddle(blocksAt, &above, &toCheck);
			while (doLoop) {
				reader = toCheck.getReader(entrySize);
				differ = toFind->compare(reader);
				if (differ < 0) {
					//  have a new lower bound for above
					above.setLocation(&toCheck);
					doLoop = below->getMiddle(blocksAt, &above, &toCheck);
				} else if (differ == 0) {
					//  match ends work
					toFind->found = true;
					doLoop = false;
				} else {
					// have better bound for below
					below->setLocation(&toCheck);
					doLoop = below->getMiddle(blocksAt, &above, &toCheck);
				}
			}
			if (!(toFind->found)) {
				// below will not need to be checked by
				//    any more targets
				// above is set will smallest range location to check
				below->setLocation(&above);
			}
		}
	}
	return (toFind->found);
}


// called with largest desired target
//    requires entries in this Blocks object to be in ascending order
//    the entry at belowSet is known to be smaller than entry
//      described by toFind
// if toFind is not found its location will describe the location
//   of the largest Blocks entry that is smaller than toFind
bool Blocks::chopTop(LinkFinder *belowSet, LinkFinder* toFind) {
	toFind->found = false;
	char *smallest = belowSet->getReader(entrySize);
	int differ = toFind->compare(smallest);
	if (differ < 0) {
		// entire range is larger than to find
		//   so set flags to end search
		toFind->belowAll = true;
		toFind->aboveAll = false;
	} else if (differ == 0) {
		// no more work if match is found
		toFind->found = true;
	} else {
		//  smallest is less than toFind
		toFind->belowAll = false;
		// create room to mark bounds for binary search
		LinkFinder above;
		LinkFinder below;
		// assume binary search will not be required
		bool doLoop = false;
		if (lastBlock == 0) {
			// largest entry is in same block as belowSet
			char *reader = belowSet->blockStart;
			reader += lastTop - entrySize;
			differ = toFind->compare(reader);
			if (0 < differ) {
				// largest smaller than toFind
				//    so it describes the desired return location
				//   and toFind will not be found
				toFind->aboveAll = true;
				toFind->belowAll = false;
				toFind->block = 0;
				toFind->blockStart = smallest;
				unsigned entryCount = lastTop / entrySize;
				toFind->inBlock = entryCount - 1;
				toFind->entriesInBlock = entryCount;
			} else if (differ == 0) {
				// no more work if match is found
				toFind->found = true;
			} else {
				// need bounds for binary search
				unsigned entryCount = lastTop / entrySize;
				below.setLocation(belowSet);
				above.setLocation(0, entryCount - 1, smallest, entryCount);
				doLoop = true;
			}
		} else {
			// get access to largest
			char *biggest = (char *)(blocksAt[lastBlock]);
			char *reader = biggest;
			reader += lastTop - entrySize;
			differ = toFind->compare(reader);
			if (0 < differ) {
				// largest smaller than toFind
				//    so it describes the desired return location
				//   and toFind will not be found
				toFind->aboveAll = true;
				toFind->belowAll = false;
				toFind->block = lastBlock;
				toFind->blockStart = biggest;
				unsigned entryCount = lastTop / entrySize;
				toFind->inBlock = entryCount - 1;
				toFind->entriesInBlock = entryCount;
			} else if (differ == 0) {
				// no more work if match is found
				toFind->found = true;
			} else {
				// need bounds for binary search
				unsigned entryCount = lastTop / entrySize;
				below.setLocation(belowSet);
				above.setLocation(
						lastBlock, entryCount - 1, biggest, entryCount);
				doLoop = true;
			}
		}
		if (doLoop) {
			// entry at below is smaller than toFind
			toFind->aboveAll = false;
			//  get middle to aid search an bounds reduction
			LinkFinder toCheck;
			doLoop = below.getMiddle(blocksAt, &above, &toCheck);
			while (doLoop) {
				char *reader = toCheck.getReader(entrySize);
				differ = toFind->compare(reader);
				if (differ < 0) {
					// reduce bounds on upper side
					above.setLocation(&toCheck);
					doLoop = below.getMiddle(blocksAt, &above, &toCheck);
				} else if (differ == 0) {
					// found match so work is over
					toFind->found = true;
					doLoop = false;
				} else {
					// reduce bounds on lower side
					below.setLocation(&toCheck);
					doLoop = below.getMiddle(blocksAt, &above, &toCheck);
				}
			}
			if (!(toFind->found)) {
				toFind->setLocation(&below);
			}
		}
	}
	return (toFind->found);
}
//  toFind describes the largest desired target that has not been
//     eliminated by a search failure
//  the belowModel to above range will be checked to see if it contains
//     an entry that matches toFind.
//  If no matching entry exists the potential exist to reduce this
//     range for future searches by adjusting the location of above
//  true is returned if a match for toFind is found
bool Blocks::chopUpper(LinkFinder* toFind, LinkFinder *belowModel,
		LinkFinder *above) {
	LinkFinder below;
	below.setLocation(belowModel);

	char *reader = above->getReader(entrySize);
	int differ = toFind->compare(reader);
	toFind->found = (differ == 0);
	toFind->belowAll = false;
	// if (differ <= 0) then either the match was found
	//    or toFind could not reduce the block searching range
	if (differ < 0) {
		// differ will cut bounds
		toFind->aboveAll = false;
		reader = below.getReader(entrySize);
		differ = toFind->compare(reader);
		if (differ < 0) {
			// the remaining range is larger than toFind
			//    so the searching can terminate in failure
			toFind->belowAll = true;
		} else if (differ == 0) {
			//  search terminated in success
			toFind->found = true;
		} else {
			// below is a possible return
			LinkFinder toCheck;
			//  chop remaining range in half to search
			//    or a large valid setting for below
			bool doLoop = below.getMiddle(blocksAt, above, &toCheck);
			while (doLoop) {
				char *reader = toCheck.getReader(entrySize);
				differ = toFind->compare(reader);
				if (differ < 0) {
					// toCheck is not a valid new setting for below
					//    but in does reduce the range for
					//    for finding the best value
					above->setLocation(&toCheck);
					doLoop = below.getMiddle(blocksAt, above, &toCheck);
				} else if (differ == 0) {
					//  a match was found so it is time to exit
					toFind->found = true;
					doLoop = false;
				} else {
					// place a larger value in below for the return
					below.setLocation(&toCheck);
					doLoop = below.getMiddle(blocksAt, above, &toCheck);
				}
			}
			if (!(toFind->found)) {
				above->setLocation(&below);
			}
		}
	}
	return (toFind->found);
}

void Blocks::absorb(Blocks *discard) {
	// make sure there is something to absorb
	if (0 < (discard->lastTop)) {
		// save address of last block for this object
		void *thisLast = blocksAt[lastBlock];

		// compute minimum # of blocks required to combine the 2 objects
		unsigned roomCount = 1 + lastBlock + discard->lastBlock;
		//  check to see if 2 final blocks can be combined in a single
		//   block
		bool overFlow = (blockSize < (lastTop + discard->lastTop));
		if (overFlow) {
			// the 2 final blocks must be store in 2 blocks
			++roomCount;
		}
		//  make sure current block list has enough places
		//    to describe all of the combined blocks
		if (atRoom < roomCount) {
			blocksAt = (void **)realloc(blocksAt, roomCount * sizeof(void *));
			atRoom = roomCount;
		}

		// get access to last block of discard
		char *otherLast = (char *)((discard->blocksAt)[discard->lastBlock]);
		if (0 < discard->lastBlock) {
			// place all full blocks from discard in the block list for this
			memcpy(blocksAt + lastBlock, discard->blocksAt,
					sizeof(void *) * (discard -> lastBlock));
			lastBlock += discard->lastBlock;
			// discard is going to be emptied but it still
			//  needs malloced memory in its first block
			//  the contents of otherLast will be moved to
			//    to another block in this object but its
			//    memory allocation will be left with discard
			(discard->blocksAt)[0] = otherLast;

			// re-insert initial tail block for this at new place in list(
			blocksAt[lastBlock] = thisLast;
		}
		// code now needs to combine data in last blocks
		//   need to add to block currently ast lastBlock
		char *tailAdd = (char *)thisLast;
		//  save data from this object that is in this block
		tailAdd += lastTop;
		if (overFlow) {
			// compute the amount of data from otherLast
			//   that can be placed in the block currently at lastBlock
			unsigned fillRoom = blockSize - lastTop;

			// prep lastTop with amount in discard's last block
			lastTop = discard->lastTop;
			if (0 < fillRoom) {
				// complete fill of current last block
				memcpy(tailAdd, otherLast, fillRoom);
				// advance otherLast to get untransfered data
				otherLast += fillRoom;
				// set lastTop with the amount of data to absorb
				//     after otherLast
				lastTop -= fillRoom;
			}
			//  a new block malloc is required to contain the overflow
			tailAdd = (char *)(malloc(blockSize));
			//  move overflow into block that will be used for
			//    new last block in this object
			memcpy(tailAdd, otherLast, lastTop);

			// set new last block
			++lastBlock;
			blocksAt[lastBlock] = (void *)tailAdd;
		} else {
			// all data in last blocks of both objects can
			//   be combined into a single block
			memcpy(tailAdd, otherLast, discard->lastTop);
			lastTop += discard->lastTop;
		}
		// indicate that discard is empty
		discard->lastBlock = 0;
		discard->lastTop = 0;
	}
}


Blocks::~Blocks() {
	for(unsigned i = 0; i <= lastBlock; i++) {
		free(*(blocksAt + i));
	}
	free(blocksAt);
}
