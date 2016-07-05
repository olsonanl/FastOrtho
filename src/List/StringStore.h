/*
 * StringStore.h
 *
 *  Created on: May 14, 2010
 *      Author: mscott
 */

#ifndef STRINGSTORE_H_
#define STRINGSTORE_H_

class StringStore {
private:
	unsigned blockSize;
	unsigned lastBlock;
	unsigned lastTop;
	void **blocksAt;
	unsigned atRoom;

public:
	StringStore(unsigned size);
	char* add(const char *toStore);
	char* debug(const char *toStore);
	char* add(const char *start, const char *stop);
	void releaseTail();
	virtual ~StringStore();
};

#endif /* STRINGSTORE_H_ */
