/*
 * Ortho.cpp
 *
 *  Created on: Jun 21, 2010
 *      Author: mscott
 */

#include "Ortho.h"
#include "Classify.h"
#include <float.h>
#include <iostream>

Ortho::Ortho(LinkParser *parser, Links** bins) {
	// get module access for link storing
	savers = bins;

	// create values to support quirky intra e value selection
	intraSave = new Blocks(sizeof(unsigned) + sizeof(float));
	intraCount = 0;

	// need array for finding best links between a query
	//   and a gene in each taxon
	taxonCount = (parser->ggParser)->getTaxonCount();
	minInters = (float *)malloc(taxonCount * sizeof(float));

	// used to block intra taxon links that are weaker than
	//   some inter taxon link

	// create array to prevent storing multiple links
	//   between a query and subject pair
	bitBytes = (parser->ggParser)->getGeneCount();
	bitBytes += 7;
	bitBytes >>= 3;
	saved = (unsigned char *)malloc(bitBytes * sizeof(unsigned char));

	unsigned loopCount = 0;

	// loop through link lines
	while((parser->rightGene) != NULL) {
		++loopCount;
		if (parser->atNewQuery) {
			// store intra links selected using quirky logic
			storeIntra();
				// deal with initial link from next query block
			setNewQuery(parser);
		}
		float nextE = parser->eValue;
		subject = (parser->rightGene)->order;
		subjectTaxon = (parser->rightGene)->parent;

		if (queryTaxon == subjectTaxon) {
			if (!(parser->isValid())) {
				if (parser->hasBadPi()) {
					// part of quirky behavior is that
					//    an invalid intra link blocks the
					//    remainder of the intra links in the query block
					minInter = nextE - 1.0f;
				}
			} else if (nextE <= minInter) {
				// save intra link information
				//    for later quirky selection
				char *storeAt = intraSave->getAdd();
				*((unsigned *)storeAt) = subject;
				storeAt += sizeof(unsigned);
				*((float *)storeAt) = nextE;
				++intraCount;
			}
		} else if (parser->isValid()){
			// check to see if inter link is a repeat
			unsigned char bit = 1;
			bit <<= (subject & 7);
			unsigned index = subject >> 3;
			if ((bit & saved[index]) == 0) {
				// first time in query block at subject
				if (minInters[subjectTaxon] < nextE) {
					addLine(OTHER_INTER_AT, nextE);
				} else {
					// this is among the strongest links
					//   to this taxon
					addLine(BEST_INTER_AT, nextE);
					// update linking level to taxon
					minInters[subjectTaxon] = nextE;
					if (nextE < minInter) {
						// update intra link blocking
						minInter = nextE;
					}
				}
				// block anymore query block links to subjectTaxon
				saved[index] |= bit;
			}
		} else {
			if (nextE < minInters[subjectTaxon]) {
				minInters[subjectTaxon] = nextE;
				if (nextE < minInter) {
					// update intra link blocking
					minInter = nextE;
				}
			}
		}
		parser->setNextLine();
	}
	// save any delayed intra links
	storeIntra();
}


int compareOrtho(const char *left, const char *right) {
	unsigned leftAt = *((unsigned *)left);
	unsigned rightAt = *((unsigned *)right);
	int result = 0;
	if (leftAt < rightAt) {
		--result;
	} else if (rightAt < leftAt) {
		++result;
	} else {
		left += sizeof(unsigned);
		float leftE = *((float*)left);
		right += sizeof(unsigned);
		float rightE = *((float *)right);
		// do quirky sorting based on e value
		if (leftE < rightE) {
			++result;
		} else if (rightE < leftE) {
			--result;
		}
	}
	return result;
}

void Ortho::storeIntra() {
	if (0 < intraCount) {
		if (1 < intraCount) {
			intraSave->sort(compareOrtho);
		}
		// march through list
		BlockMarch * marcher = intraSave->getFullMarch();
		char *reader = marcher->getEntry();
		// always save first entry
		subject = *((unsigned *)reader);
		reader += sizeof(unsigned);
		addLine(TOP_INTRA_AT, *((float *)reader));
		reader = marcher->getNextEntry();
		while(reader != NULL) {
			// keep links unique by blocking
			//   repeated subjects
			unsigned nextSubject = *((unsigned *)reader);
			if (nextSubject != subject) {
				subject = nextSubject;
				reader += sizeof(unsigned);
				addLine(TOP_INTRA_AT, *((float *)reader));
			}
			reader = marcher->getNextEntry();
		}
		// empty ram for next query block
		intraSave->empty();
		intraCount = 0;
	}
}

void Ortho::addLine(unsigned type, float eValue) {
	if (query < subject) {
		(savers[type])->add(query, subject, eValue);
	} else {
		(savers[type])->add(subject, query, eValue);
	}
}

void Ortho::setNewQuery(LinkParser *parser) {
	// reset to taxon blocks to for new query
	std::fill(minInters, minInters + taxonCount, FLT_MAX);
	std::fill(saved, saved + bitBytes, (unsigned char)0);

	// collect information about link
	query = (parser->leftGene)->order;
	queryTaxon = (parser->leftGene)->parent;
	subject = (parser->rightGene)->order;
	subjectTaxon = (parser->rightGene)->parent;
	minInter = FLT_MAX;
}


Ortho::~Ortho() {
	free(minInters);
	free(saved);
	delete(intraSave);
}
