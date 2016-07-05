/*
 * LinePull.cpp
 *
 *  Created on: Feb 26, 2010
 *      Author: mscott
 */

#include "LinePull.h"
#include <string>
#include <vector>
#include <ctype.h>
#include <iostream>
#include <stdio.h>
#include <string.h>


LinePull::LinePull(size_t room, const char *filePath) {
	pullTries = 0;
	bufferRoom = room;
	partRoom = 0;
	reader = fopen64(filePath, "r");
	if (reader == NULL) {
		// block any attempt to free buffers
		partPoints = NULL;
		buffer = NULL;
		// indicate that nothing is available
		available = 0;
	} else {
		// create buffer for loading data form file
		buffer = (char *)malloc(bufferRoom);
		// mark end of buffer
		pastBuffer = buffer + bufferRoom;
		// set virgin data pointer ast start of buffer
		unused = buffer;
		// attempt to fill buffer
		available = fread(buffer, sizeof(char), bufferRoom, reader);
		if (available < bufferRoom) {
			// string aid when buffer is not completely fill
			buffer[available] = '\0';
		}
		//  create pointer array for parts in line split work
		partRoom = 0x100;
		partPoints = (char **)(malloc(partRoom * sizeof(char *)));
	}
}

bool LinePull::refreshBuffer() {
	// need something to discard or buffer can not be refreshed
	bool result = (buffer < unused);
	if (result) {
		// if initial load or last refresh must have filled buffer
		//    or there is no more file data for refresh
		result = (pastBuffer <= unused + available);
		// shift to free room for load
		memmove(buffer, unused, available);
		unused = buffer;
		if (result) {
			// try to load more data
			size_t readCount =
				fread(buffer + available, sizeof(char), bufferRoom - available,
					 reader);
			// result is true only if more data is obtained
			result = (0 < readCount);
			available += readCount;
		}
	}
	return result;
}

char ** LinePull::getCharSplitParts(char splitter, unsigned & partCount) {
	++pullTries;
	char **result = NULL;
	partCount = 0;
	// skip white space at start of line
	while ((0 < available) && (*unused != splitter) && (isspace(*unused))) {
		++unused;
		--available;
		if (available == 0) {
			refreshBuffer();
		}
	}
	// make sure full line will fit into buffer
	char *lineEnd = (char *)memchr(unused, '\n', available);
	if (lineEnd == NULL) {
		if (refreshBuffer()) {
			lineEnd = (char *)memchr(unused, '\n', available);
		}
	}
	if (lineEnd == NULL) {
		// could not find termination so assume termination
		//  at end of file data
		lineEnd = unused + available;
		if (pastBuffer <= lineEnd) {
			// did not reach end of file data so terminate file work
			available = 0;
		}
	}
	if (0 < available) {
		while (unused <= lineEnd) {
			// find end of next part
			char *splitAt = (char *)memchr(unused, splitter, lineEnd - unused);
			if (splitAt == NULL) {
				// splitting character was not available so
				//   start work with line as a single part
				splitAt = lineEnd;
			}
			// trim white front of part
			while((unused < splitAt) && isspace(*unused)) {
				--available;
				++unused;
			}
			// make sure partPoints has enough room
			//   to hold another part
			if (partRoom <= partCount) {
				partRoom += 0x100;
				partPoints =
					(char **)(realloc(partPoints, partRoom * sizeof(char *)));
			}
			// store next in partPoints
			char *start = unused;
			partPoints[partCount] = start;
			++partCount;

			// buffer used to hold part is now considered used
			unused = splitAt + 1;
			if (available < unused - start) {
				available = 0;
			} else {
				available -= unused - start;
			}

			// trim white tail of part
			--splitAt;
			while ((start <= splitAt) && (isspace(*splitAt))) {
				--splitAt;
			}
			// splitAt is either less than start or pointing at non-white
			++splitAt;
			// move splitAt back to whitepace and string terminate part
			*splitAt = '\0';
		}
		result = partPoints;
	}
	return result;
}

char* LinePull::trim(char *start, char *stop) {
	while ((isspace(*start)) && (start < stop)) {
		++start;
	}
	while ((start < stop) && (isspace(*(stop - 1)))) {
		--stop;
	}
	*stop = '\0';
	return start;
}

char* LinePull::trim(char *start) {
	return (trim(start, start + strlen(start)));
}


char * LinePull::getFrontWhiteSplit(){
	char *result = NULL;
	// search for whitespace skipping empty lines
	while ((0 < available) && (isspace(*unused))) {
		++unused;
		--available;
		if (available == 0) {
			refreshBuffer();
		}
	}
	if (0 < available) {
		if (available == 1) {
			refreshBuffer();
		}
		// found some non-white space
		char *whiteFind = unused + 1;
		char *reload = unused + available;
		// search for end of white space
		while ((whiteFind < reload) && (!isspace(*whiteFind))) {
			++whiteFind;
			//  adjust if buffer is empty and search did not examine
			// an entire buffer's work of data
			if ((whiteFind == reload) && (buffer < unused)) {
				// adjust whiteFind for shift that will be
				//    generated by buffer refresh
				whiteFind -= unused - buffer;
				refreshBuffer();
				reload = unused + available;
			}
		}
		if (whiteFind < reload) {
			// found a whitespace character so check for end of line
			didLastWhite = (*whiteFind == '\n');
			// string terminate non-white block from unused to whiteFind
			*whiteFind = '\0';
			// new unused will be just past the termination
			++whiteFind;
			if (reload <= whiteFind) {
				// termination was at end of buffer
				//    a refresh is required to prevent available
				//    going to 0 thus indicating end of file work

				// adjust whiteFind for change after refresh shift
				whiteFind -= unused - buffer;
				refreshBuffer();
				// set new reload
				reload = unused + available;
			}
			// set return at start of non-white block
			result = unused;
			// remove non-white string from unused portion of buffer
			unused = whiteFind;
			//  adjust available for new unused pointer
			available = reload - unused;
		} else {
			// no non-whitespace in remaining file data
			if (whiteFind < pastBuffer) {
				// ran out of file data
				*whiteFind = '\0';
				result = unused;
			} //else {
				// entire buffer is full on non whitespace characters
				//   treat remainder of file as invalid and
				//   leave result as NULL for invalid return
			//}
			// block anymore work on file
			available = 0;
		}
	}
	return result;
}

char * LinePull::getNextWhiteSplit(){
	char *result = NULL;
	// discard whitespace at from of line remainder at unused
	while ((!didLastWhite) && (0 < available) && (isspace(*unused))) {
		// check for end of line
		didLastWhite = (*unused == '\n');
		// skip over whitespace
		++unused;
		// adjust to match skip
		--available;
		if (available == 0) {
			refreshBuffer();
		}
	}
	// check to see if there was any non-whitespace in line
	//   remainder stored in buffer
	if ((!didLastWhite) && (0 < available)) {
		if (available == 1) {
			refreshBuffer();
		}
		// non-whitespace at unused so start white search just after it
		// start search for whitespace
		char *whiteFind = unused + 1;
		// aid to keep whiteFind in actual data
		char *reload = unused + available;
		// search for end of non-whitespace
		//    front is at unused
		while ((whiteFind < reload) && (!isspace(*whiteFind))) {
			++whiteFind;
			//  adjust if buffer is empty
			//     and search did not examine and entire buffer's work
			//     of data
			if ((whiteFind == reload) && (buffer < unused)) {
				// adjust whiteFind for shift that will be
				//    generated by buffer refresh
				whiteFind -= unused - buffer;
				refreshBuffer();
				reload = unused + available;
			}
		}
		if (whiteFind < reload) {
			// found whitespace to terminate return
			//   if this is last non-white part in line set the flast
			didLastWhite = (*whiteFind == '\n');
			// terminate non-white block starting at unused
			*whiteFind = '\0';
			// prepare to remove this non-white block from unused data
			++whiteFind;
			if (reload <= whiteFind) {
				// termination at end of buffer can lead to premature view
				//   that file data is gone, so a refresh is in order
				whiteFind -= unused - buffer;
				refreshBuffer();
				reload = unused + available;
			}
			result = unused;
			available -= whiteFind - result;
			unused = whiteFind;
		} else {
			// no non-whitespace in remaining file data
			if (whiteFind < pastBuffer) {
				// ran out of file data
				*whiteFind = '\0';
				result = unused;
			} //else {
				// entire buffer is full on non whitespace characters
				//   treat remainder of file as invalid and
				//   leave result as NULL for invalid return
			//}
			// block anymore work on file
			available = 0;
		}
	}
	return result;
}



char *LinePull::getLine() {
	char *result = NULL;
	// never exit with available == 0 unless no hope for next getLine()
	if (0 < available) {
		char *lineEnd = (char *)memchr(unused, '\n', available);
		if (lineEnd == NULL) {
			// see buffer refresh will provide access to line end
			refreshBuffer();
			lineEnd = (char *)memchr(unused, '\n', available);
		}
		if (lineEnd != NULL) {
			// found line end
			if  (pastBuffer <= lineEnd + 1) {
				// line end is at end of buffer
				//  block confusion by refreshing buffer
				lineEnd -= unused - buffer;
				refreshBuffer();
			}
			result = unused;
			// mark data at end of line
			unused = lineEnd + 1;
			available -= unused - result;
		} else {
			// no '\n' in full buffer or remainder of file
			lineEnd = unused + available;
			if (lineEnd < pastBuffer) {
				result = unused;
			} else {
				--lineEnd;
			}
			// this is the last line
			available = 0;
		}
		// replace line termination with zeros
		*lineEnd = '\0';
		--lineEnd;
		// discard Microsoft line termination
		if ((result <= lineEnd) && (*lineEnd == '\r')) {
			*lineEnd = '\0';
		}
	}
    return result;
}

LinePull::~LinePull() {
	if (reader != NULL) {
		free(buffer);
		free(partPoints);
		fclose(reader);
		reader = NULL;
	}
}
