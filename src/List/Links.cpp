/*
 * Links.cpp
 *
 *  Created on: Apr 9, 2010
 *      Author: mscott
 */

#include "Links.h"
#include <iostream>
#include <stdio.h>
#include <math.h>
#include <string.h>

Links::Links() {
	entrySize = sizeof(float) + 2 * sizeof(unsigned);
	storage = new Blocks(entrySize);
	intraDepths = NULL;
	intraValues = NULL;
}

Blocks* Links::getStorage() {
	return storage;
}

unsigned Links::getEntryCount() {
	return (storage->getEntryCount());
}

BlockMarch *Links::getFullMarch() {
	return (storage->getFullMarch());
}

void Links::add(unsigned left, unsigned right, float rating) {
	char* toFill = storage->getAdd();
	*((unsigned *)toFill) = left;
	toFill += sizeof(unsigned);
	*((unsigned *)toFill) = right;
	toFill += sizeof(unsigned);
	*((float *)toFill) = rating;
}

void Links::add(char *toAdd) {
	storage->add(toAdd);
}

void Links::absorb(Links* discard) {
	storage->absorb(discard->storage);
}

// this should only be called after adjustIntras which is
//   responsible for setting taxonCount
// and creating intraDepths, and intraValues
void Links::moveIntras(Links* bestInters, char *linkTrack) {
	for (unsigned i = 0; i < taxonCount; i++) {
		unsigned depth = intraDepths[i];
		if (0 < depth) {
			//  have intra links for taxon i
			char *reader = intraValues[i];
			// initialize aids for computing mcl weight for links
			double adjust = 0.0;
			double bestSum = 0.0;
			unsigned bestCount = 0;
			for (unsigned j = 0; j < depth; j++) {
				unsigned left = *((unsigned *)reader);
				reader += sizeof(unsigned);
				unsigned right = *((unsigned *)reader);
				reader += sizeof(unsigned);
				float value = *((float *)reader);
				reader += sizeof(float);
				adjust += value;
				if  ((linkTrack[left] != '\0')
							||
					 (linkTrack[right] != '\0')) {
					// track e-value sum for genes that
					//    link to other taxon
					bestSum += value;
					++bestCount;
				}
			}
			if (bestCount == 0) {
				// none of the genes linked to other taxon
				//   we will adjust with the mean of all e-values
				//   for genes in taxon i
				adjust /= depth;
			} else {
				//  use mean of weights for links to genes
				//    that also link outside of taxon i
				adjust = bestSum / bestCount;
			}
			// reset to allow march through all of intra links for taxon i
			reader = intraValues[i];
			while (0 < depth) {
				// gene indexes will not change
				unsigned left = *((unsigned *)reader);
				reader += sizeof(unsigned);
				unsigned right = *((unsigned *)reader);
				reader += sizeof(unsigned);
				float value = *((float *)reader);
				reader += sizeof(float);
				// e-value will be divided by a mean
				//  all e-values will be negative so the division
				//    will result in the desired sign change
				value = (float)( value /adjust);
				// save with bestInters to form final single list
				//   for output
				bestInters->add(left, right, value);
				--depth;
			}
			// release ram that is no longer required
			free(intraValues[i]);
		}
	}
	free(intraDepths);
	intraDepths = NULL;
	free(intraValues);
	intraValues = NULL;
}

// standard sorting based on first coordinate being most significant
int compare(const char *left, const char *right) {
	int result = 0;
	unsigned seeLeft = ((unsigned *)left)[0];
	unsigned seeRight = ((unsigned *)right)[0];
	if (seeLeft < seeRight) {
		--result;
	} else if (seeRight < seeLeft) {
		++result;
	} else {
		seeLeft = ((unsigned *)left)[1];
		seeRight = ((unsigned *)right)[1];
		if (seeLeft < seeRight) {
			--result;
		} else if (seeRight < seeLeft) {
			++result;
		}
	}
	return result;
}

void Links::sort() {
	storage->sort(compare);
}

bool merge(const char *left, const char *right, char *saveIn) {
	// check to see if both coordinates of left & right match
	bool result = (0 == memcmp(left, right, 2 * sizeof(unsigned)));
	if (result) {
		// coordinates do not change
		memcpy(saveIn, left, 2 * sizeof(unsigned));
		// adjust to work with e-values
		left += 2 * sizeof(unsigned);
		right += 2 * sizeof(unsigned);
		saveIn += 2 * sizeof(unsigned);
		// save mean of e-values
		float value = *((float *)left);
		value += *((float *)right);
		value *= 0.5f;
		*((float *)saveIn) = value;
	}
	return result;
}

void Links::keepPairs() {
	storage->sort(compare);
	storage->keepPairs(merge);
}

void Links::demoteSingles(Links *lesser, char *linkTrack) {
	storage->sort(compare);
	storage->demoteSingles(merge, lesser->storage, linkTrack);
}

void Links::adjustForIntra(unsigned genomeCount, unsigned* toParent) {
	getParent = toParent;
	taxonCount = genomeCount;
	// prepare array that will count number of intra links per taxon
	intraDepths = (unsigned *)(malloc(taxonCount * sizeof(unsigned)));
	// prepare array that will store intra links per taxon
	intraValues = (char **)(malloc(taxonCount * sizeof(char *)));
	for (unsigned i = 0; i < taxonCount; i++) {
		intraDepths[i] = 0;
	}

	// aid for collecting per taxon data
	Blocks collector(entrySize);

	// set atParent to non-existent taxon to trigger new taxon
	//   response to first taxon discovered
	unsigned atParent = taxonCount;

	// used to tra
	unsigned atCount = 0;

	//  work from top of storage so ram
	//    can be released after data is moved
	char *toDrop = storage->getLastEntry();
	while (toDrop != NULL) {
		char *endEntry = toDrop;
		unsigned low = *((unsigned *)endEntry);
		unsigned nextParent = getParent[low];
		// check to see if intra entries have transitioned to a
		//   new taxon
		if (nextParent != atParent) {
			if (0 < atCount) {
				collector.sort(compare);
				// move collector contents into a single block
				intraValues[atParent] = collector.emptyToOneBlock();
				// entries are twice atCount because of bi-directional save
				intraDepths[atParent] = (atCount << 1);
			}
			atParent = nextParent;
			atCount = 0;
		}
		++atCount;
		// add next member to collector
		collector.add(endEntry);

		//  add next member with left & right reversed
		char *storeIn = collector.getAdd();
		char *reader = endEntry;
		// advance to store right as left
		reader += sizeof(unsigned);
		memcpy(storeIn, reader, sizeof(unsigned));
		// advance storeIn and store left as right
		storeIn += sizeof(unsigned);
		memcpy(storeIn, endEntry, sizeof(unsigned));
		// advance reader to e-value then advance storeIn to save
		reader += sizeof(unsigned);
		storeIn += sizeof(unsigned);
		memcpy(storeIn, reader, sizeof(float));

		// remove toDrop from intra
		toDrop = storage->dropLastEntry(toDrop);
	}
	// dump final taxon entry block if found
	if (0 < atCount) {
		collector.sort(compare);
		intraValues[atParent] = collector.emptyToOneBlock();
		intraDepths[atParent] = (atCount << 1);
	}
}


void Links::getClones(unsigned lowMatch, std::vector<unsigned> *result) {
	unsigned parentAt = getParent[lowMatch];
	unsigned aboveIn = intraDepths[parentAt];
	bool saveMatch = true;
	// need some entries for the parentAt taxon
	if (0 < aboveIn) {
		//  get entries for parentAt taxon
		char *entries = intraValues[parentAt];
		// check to see if 1st exceeds any lowMatch entries
		//   which will block discovery of any entries
		unsigned bottom = *((unsigned *)entries);
		if (bottom <= lowMatch) {
			//  range of entries might contain lowMatch
			//    next objective is to set startIn with smallest
			//      entry that equals lowMatch
			unsigned startIn = 0;
			// will also use insideIn as starting point for
			//    finding end of entries which match lowMatch
			unsigned insideIn = 0;
			//  if bottom matches then startIn is known
			//      and insideIn is a valid value
			if (bottom < lowMatch) {
				//  set insideIn to indicate that no matches have been found
				insideIn = aboveIn;
				//  set belowIn with location that is too small to match lowMatch
				unsigned belowIn = 0;
				//  check midway between known invalid locations
				unsigned checkIn = aboveIn >> 1;
				while (belowIn < checkIn) {
					//  get value at
					char *reader = entries + (checkIn * entrySize);
					unsigned low = *((unsigned *)reader);
					if (low < lowMatch) {
						// new value for lower failures
						belowIn = checkIn;
						checkIn = (belowIn + aboveIn) >> 1;
					} else if (low == lowMatch) {
						// save location of discovered match
						insideIn = checkIn;
						//  set checkIn to exit loop
						checkIn = belowIn;
					} else {
						// new location for upper failures
						aboveIn = checkIn;
						checkIn = (belowIn + aboveIn) >> 1;;
					}
				}
				// check to see lowMatch was found
				if (insideIn < aboveIn) {
					//  next set startIn with smallest index that points
					//    at a lowMatch entry
					// start with a known valid entry
					startIn = insideIn;
					checkIn = (belowIn + startIn) >> 1;
					while (belowIn < checkIn) {
						char *reader = entries + (checkIn * entrySize);
						unsigned low = *((unsigned *)reader);
						if (low < lowMatch) {
							//  new location for lower failure
							belowIn = checkIn;
						} else {
							// move startIn to valid smaller value
							startIn = checkIn;
						}
						// move to middle
						checkIn = (belowIn + startIn) >> 1;
					}
				}
			}
			if (insideIn < aboveIn) {
				//  startIn is properly set
				//    need to move aboveIn to smallest non-matching index
				unsigned checkIn = (insideIn + aboveIn) >> 1;
				while (insideIn < checkIn) {
					char *reader = entries + (checkIn * entrySize);
					unsigned low = *((unsigned *)reader);
					if (lowMatch < low) {
						// decrease upper bound for final aboveIn
						aboveIn = checkIn;
					} else {
						// increase lower bound for final aboveIn
						insideIn = checkIn;
					}
					checkIn = (insideIn + aboveIn) >> 1;
				}
				//  return high value for all entries that match lowMatch
				char *reader = entries + (startIn * entrySize);
				// skip over portion which matches lowMatch
				reader += sizeof(unsigned);
				do {
					// get value paired with lowMatch
					unsigned next = *((unsigned *)reader);
					// check for proper place to insert lowMatch
					if ((saveMatch) && (lowMatch < next)) {
						// at first range value larger than lowMatch
						result->push_back(lowMatch);
						saveMatch = false;
					}
					result->push_back(next);
					reader += entrySize;
					++startIn;
				} while (startIn < aboveIn);
			}
		}
	}
	if (saveMatch) {
		result->push_back(lowMatch);
	}
}


Links::~Links() {
	if (intraDepths != NULL) {
		free(intraDepths);
	}
	if (intraValues != NULL) {
		free(intraValues);
	}
	delete (storage);
}
