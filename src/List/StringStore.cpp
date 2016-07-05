/*
 * StringStore.cpp
 *
 *  Created on: May 14, 2010
 *      Author: mscott
 */

#include "StringStore.h"
#include <stdlib.h>
#include <string.h>
#include <iostream>
#include <stdio.h>

StringStore::StringStore(unsigned size) {
	lastBlock = 0;
	lastTop = 0;
	blockSize = size;
	atRoom = 0x100;
	blocksAt = (void **)malloc(atRoom * sizeof(void *));
	*blocksAt = malloc(blockSize);
}

char* StringStore::add(const char *toStore) {
	unsigned length = strlen(toStore);
	char *result;
	//check to see if current block can hold copy of toStore
	if (lastTop + length < blockSize) {
		// set to store in current block
		result = (char *)blocksAt[lastBlock];
		result += lastTop;
	} else {
		// need a new block
		++lastBlock;
		// make sure new block can be tracked
		if (atRoom <= lastBlock) {
			atRoom += 0x100;
			blocksAt = (void **)realloc(blocksAt, atRoom * sizeof(void *));
		}
		// save location of new block
		blocksAt[lastBlock] = malloc(blockSize);
		// set to store at start of new block
		result = (char *)(blocksAt[lastBlock]);
		lastTop = 0;
	}
	// allow for storage of terminating '\0'
	++length;
	memcpy(result, toStore, length);
	// set lastTop for next storage location in current block
	lastTop += length;
	return result;
}


char *StringStore::add(const char *start, const char *stop) {
	unsigned length = stop - start;
	char *result;
	if (lastTop + length < blockSize) {
		// set result to store in current block
		result = (char *)blocksAt[lastBlock];
		result += lastTop;
	} else {
		// need a new block
		++lastBlock;
		if (atRoom <= lastBlock) {
			// need more room for storing block locations
			atRoom += 0x100;
			blocksAt = (void **)realloc(blocksAt, atRoom * sizeof(void *));
		}
		// get ram for next block
		blocksAt[lastBlock] = malloc(blockSize);
		result = (char *)(blocksAt[lastBlock]);
		lastTop = 0;
	}
	// store desired text
	memcpy(result, start, length);
	// add termination to create a string
	result[length] = '\0';
	// set lastTop for next storage location
	lastTop += length + 1;
	return result;
}

void StringStore::releaseTail() {
	// relase unused memory in last block
	blocksAt[lastBlock] = realloc(blocksAt[lastBlock], lastTop);
}

StringStore::~StringStore() {
	// release allocated blocks of memory
	for (unsigned i = 0; i <= lastBlock; i++) {
		free(blocksAt[i]);
	}
	// release location storage for memory blocks
	free(blocksAt);
}
