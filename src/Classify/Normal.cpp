/*
 * normal.cpp
 *
 *  Created on: Jun 21, 2010
 *      Author: mscott
 */

#include "Normal.h"
#include "Classify.h"
#include <float.h>
#include <iostream>

Normal::Normal(LinkParser *parser, Links** bins) {
	savers = bins;

	// need array for finding best links between a query
	//   and a gene in each taxon
	unsigned taxonCount = (parser->ggParser)->getTaxonCount();
	float* minInters = (float *)malloc(taxonCount * sizeof(float));

	// used to block intra taxon links that are weaker than
	//   some inter taxon link
	float minInter = FLT_MAX;

	// create array to prevent storing multiple links
	//   between a query and subject pair
	unsigned bitBytes = (parser->ggParser)->getGeneCount();
	bitBytes += 7;
	bitBytes >>= 3;
	unsigned char *saved =
		(unsigned char *)malloc(bitBytes * sizeof(unsigned char));
	while((parser->rightGene) != NULL) {
		//checkSpecial(parser);
		if (parser->atNewQuery) {
			// reset convenience ram
			query = (parser->leftGene)->order;
			queryTaxon = (parser->leftGene)->parent;
			// reset testing ram
			minInter = FLT_MAX;
			std::fill(minInters, minInters + taxonCount, FLT_MAX);
			std::fill(saved, saved + bitBytes, (unsigned char)0);
		}
		nextE = parser->eValue;
		subject = (parser->rightGene)->order;
		//  first time in block
		//    now classify link before storing
		unsigned subjectTaxon = (parser->rightGene)->parent;
		if (parser->isValid()) {
			// check to see if remainder of query block should be skipped
			// check to see if this is the first time
			//   the subject has appeared in the query block
			unsigned char bit = 1;
			bit <<= (subject & 7);
			unsigned index = subject >> 3;
			if ((bit & saved[index]) == 0) {
				if (queryTaxon == subjectTaxon) {
					//  intra links must be as strong as all
					//     inter links or storage is blocked
					if (nextE <= minInter) {
						addLine(TOP_INTRA_AT);
					}
				} else {
					// check to see if link is the
					//    best between query and subjectTaxon
					if (minInters[subjectTaxon] < nextE) {
						// if not classify as other
						addLine(OTHER_INTER_AT);
					} else {
						addLine(BEST_INTER_AT);
						//  update inter taxon block
						minInters[subjectTaxon] = nextE;
						if (nextE < minInter) {
							// update intra block
							minInter = nextE;
						}
					}
				}
				// block storage of non-unique links
				saved[index] |= bit;
			}
		} else if ((queryTaxon != subjectTaxon)
						&&
					(nextE < minInters[subjectTaxon])){
			minInters[subjectTaxon] = nextE;
		}
		parser->setNextLine();
	}
	free(minInters);
	free(saved);
}

void Normal::addLine(unsigned type) {
	// store with gene indexes in smaller, larger order
	if (query < subject) {
		(savers[type])->add(query, subject, nextE);
	} else {
		(savers[type])->add(subject, query, nextE);
	}
}

Normal::~Normal() {
}
