/*
 * OrthoTop.cpp
 *
 *  Created on: Apr 8, 2010
 *      Author: mscott
 */

#include "OrthoTop.h"
#include "Parser/GgParser.h"
#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <float.h>
#include "List/Links.h"
#include <time.h>
#include "Classify/Classify.h"
#include "Executer.h"
#include <string.h>


OrthoTop::OrthoTop(OptionPass *options) {

	// ggParse will have been set by OptionPass or
	//    areOk() funtion would block constructor call in FastOrtho
	ggParser = options->getGgParser();

	// create bins that allow links to be classified
	intraTops = new Links();
	bestInters = new Links();
	otherInters = new Links();

	// put bins in array for easy passing
	Links** savers = (Links **)malloc(SAVE_TYPE_COUNT * sizeof(Links *));
	savers[TOP_INTRA_AT] = intraTops;
	savers[BEST_INTER_AT] = bestInters;
	savers[OTHER_INTER_AT] = otherInters;

	time_t firstTime;
	time (&firstTime);

	// go read link source and put discovered links in bins
	bool ok = options->classifyLinks(savers);
	//  the array is no longer required
	free(savers);
	if (ok) {
		time_t startTime;
		time (&startTime);
		// release memory required to deal with gene names in text source
		ggParser->releaseNameRam();

		// get number of taxons for
		taxonCount = ggParser->getTaxonCount();

		// get map from gene index to genome index
		unsigned *getParent = ggParser->getParentIndex();
		unsigned geneCount = ggParser->getGeneCount();

		// create array to track genes that have optimal ortholog links
		char * bestInterLinked = (char *)(malloc(geneCount * sizeof(char)));
		for (unsigned i = 0; i < geneCount; i++) {
			bestInterLinked[i] = '\0';
		}

		// remove intra links which are not bi-directional
		intraTops->keepPairs();

		// demote best inter links which are not bi-directional
		bestInters->demoteSingles(otherInters, bestInterLinked);

		//  discard secondary ortholog links that are not bi-directional
		otherInters->keepPairs();

		// re-position paralog links to support hub work
		//  i.e. promotion of secondary ortholog links into mcl input
		intraTops->adjustForIntra(taxonCount, getParent);


		//  check to see which secondary inter-links should
		//    included in mcl input
		if (0 < intraTops->getEntryCount()) {
			expandInters();
		}

		// get values for normalizing weight of ortholog links
		setInterMeans(taxonCount, getParent);

		// normalize weights of ortholog links for mcl input
		adjustInters(getParent);

		// release getParent to conserve memory for mcl processing
		ggParser->releaseParentIndex();
		getParent = NULL;

		// combine orthlog and paralog links into a single list
		//    to use as mcl input and normalize paralog weights in the process
		intraTops->moveIntras(bestInters, bestInterLinked);

		// finished with bestInterLinked so release ram
		free(bestInterLinked);

		//  finished with intraTops so release ram
		delete(intraTops);

		bool haveData = (0 < bestInters->getEntryCount());
		if (haveData) {
			//  bestInters now contains all edges for mcl
			//    sort to simplify mcl input file preparation
			bestInters->sort();

			// create file for mcl input including the ".idx"
			//    for converting mcl output in a readable file
			writeMclInput(options->getMatrixPath());
		} else {
			printf("None of the blast hits pass OrthoMCL conditions\n");
		}

		// release ram ram
		delete(bestInters);

		time_t endTime;
		time (&endTime);
		double duration = difftime (endTime, startTime);
		printf (" %.2lf to prepare classified hits for mcl\n", duration );
		time (&startTime);
		// run mcl using prepared input file
		if ((haveData) && (options->runMcl())) {
			// if run was successful create readable version of output
			doOrthoConvert(options->getMclOutputPath(),
						   ggParser->indexFile,
						   options->resultFile);
		}
		time(&endTime);
		duration = difftime (endTime, startTime);
		printf (" %.2lf to run mcl and convert its output\n", duration );
		duration = difftime (endTime, firstTime);
		printf (" %.2lf total duration\n", duration );
	}



}


void OrthoTop::setInterMeans(unsigned taxonCount, unsigned *getParent) {
	// create and initialize arrays from computing taxon pair means
	//     dealing with inter only pairs reduce required dimension for array
	//        (lower index is not required on left)
	unsigned count = taxonCount - 1;
	interMeans = (double **)malloc(count * sizeof(double *));
	unsigned **interCounts = (unsigned **)malloc(count * sizeof(unsigned *));
	for (unsigned i = 1; i < taxonCount; i++) {
		// allocate room for row
		double *means = (double *)malloc(i * sizeof(double));
		unsigned *counts = (unsigned *)malloc(i * sizeof(unsigned));
		// initialize row values to 0
		for (unsigned j = 0; j < i; j++) {
			means[j] = 0.0;
			counts[j] = 0;
		}
		// save initialized new rows
		interMeans[i - 1] = means;
		interCounts[i - 1] = counts;
	}

	Blocks *storage = bestInters->getStorage();
	BlockMarch* marcher = storage->getFullMarch();
	char *reader = marcher->getEntry();
	while (reader != NULL) {
		unsigned low = *((unsigned *)reader);
		reader += sizeof(unsigned);
		unsigned high = *((unsigned *)reader);
		// convert gene indices into taxon indices
		low = getParent[low];
		high = getParent[high];

		//  decrease upper to match array row shift
		--high;

		//  get e value for pair
		reader += sizeof(unsigned);
		float eValue = *((float *)reader);

		// bump count form mean computation
		unsigned *counter = interCounts[high];
		++counter[low];

		// increase sum for mean computation
		double *meaner = interMeans[high];
		meaner[low] += eValue;

		reader = marcher->getNextEntry();
	}
	// release object that marched through bestInters
	delete (marcher);

	// compute means
	for (unsigned i = 1; i < taxonCount; i++) {
		unsigned *counter = interCounts[i - 1];
		double *meaner = interMeans[i - 1];
		for (unsigned j = 0; j < i; j++) {
			count = counter[j];
			if (0 < count) {
				double sum = meaner[j];
				//  compute mean for taxon pair
				sum /= count;
				meaner[j] = sum;
			}
		}
	}
	// now that means are computed ram for count values
	//   is no longer required
	for (unsigned i = 1; i < taxonCount; i++) {
		free(interCounts[i - 1]);
	}
	free(interCounts);
}

// handles promotion of secondary ortholog links to
//    inputs for mcl
void OrthoTop::expandInters() {
	// collects promoted links
	Links* keepers = new Links();

	std::vector<unsigned> lowClones;
	std::vector<unsigned> highClones;
	//start from end of secondary to ease release of ram
	//    as each entry is considered
	Blocks *toCheck = otherInters->getStorage();
	char* toDrop = toCheck->getLastEntry();
	while (toDrop != NULL) {
		char *reader = toDrop;
		unsigned low = *((unsigned *)reader);
		intraTops->getClones(low, &lowClones);
		reader += sizeof(unsigned);
		unsigned high = *((unsigned *)reader);
		intraTops->getClones(high, &highClones);
		bool keep = false;

		// a secondary will be promoted only if there exist
		//    a primary link between to genes and
		//    where the 2 genes are paralogs of the secondary ends
		if (lowClones.size() == 1) {
			// the low gene does not have any paralogs
			//   so its promotion requires a paralog for the high end
			if (1 < highClones.size()) {
				keep = doHubCheck(low, highClones);
			}
		} else {
			if (1 < highClones.size()) {
				//  the
				keep = doHubCheck(lowClones, highClones);
			} else {
				//  no paralogs for the high end but
				//    we can still use paralogs for the low end
				keep = doHubCheck(lowClones, high);
			}
		}
		if (keep) {
			keepers->add(toDrop);
		}
		// finished with toDrop so memory relase from toCheck is possible
		toDrop = toCheck->dropLastEntry(toDrop);
		// clear paralogs for next loop
		lowClones.clear();
		highClones.clear();
	}
	bestInters->absorb(keepers);
	delete(otherInters);
	delete (keepers);
}

void OrthoTop::writeMclInput(char *matrixPath) {
	unsigned geneCount = ggParser->getGeneCount();
	// bestInters contains bi directional version so of all weigthted edges
	//   where each gene is represented by the order in which it appeared
	//   in the .gg or .glg file.
	//   pairs are sorted by left gene index then right gene index

	// It may be the case that some genes have no weighted edges
	//   but it may be desirable for mcl to receive and continuous range
	//   	of indices
	// the work below will create a translate mapping to take
	//   the gene indexes with weighted edges to a continuous range
	//   starting at 0
	unsigned* translate = (unsigned *)malloc(geneCount * sizeof(unsigned));
	unsigned foundCount = 0;
	//  set lastLow to insure that smallest weighted edge is
	//  included in the mapping
	unsigned lastLow = geneCount;

	BlockMarch* marcher = bestInters->getFullMarch();
	char *reader = marcher->getEntry();

	while (reader != NULL) {
		unsigned nextLow = *((unsigned *)reader);
		if (nextLow != lastLow) {
			//  set mapping entry to take new gene to
			//   it compressed position
			translate[nextLow] = foundCount;
			// also include this gene in the .glg file for
			//   those genes with weighted edges
			ggParser->addToIndexFile(nextLow, foundCount);
			// indicate that nextLow has been mapped
			lastLow = nextLow;
			//  adjust compressed position
			++foundCount;
		}
		reader = marcher->getNextEntry();
	}
	// mapping for output is complete
	ggParser->closeIndexFile();
	if (0 < foundCount) {
		// insert standard header for input to mcl file
		FILE * writer = fopen(matrixPath, "w");
		fprintf(writer, "(mclheader\n");
		fprintf(writer, "mcltype matrix\n");
		fprintf(writer, "dimensions %dx%d\n", foundCount, foundCount);
		fprintf(writer, ")\n\n");
		fprintf(writer, "(mclmatrix\n");
		fprintf(writer, "begin\n\n");


		// loop through links generated by orthomcl filtering and logic
		marcher = bestInters->getFullMarch();
		reader = marcher->getEntry();
		//  set start of lined with next side
		unsigned lastLeft = *((unsigned *)reader);
		fprintf(writer, "%d   ", translate[lastLeft]);
		// complete mcl input line with indices linked to lastLeft
		while (reader != NULL) {
			// get next left index
			unsigned index = *((unsigned *)reader);
			if (index != lastLeft) {
				// terminate line for lastleft
				fprintf(writer, " $\n");
				// update new lastLeft
				lastLeft = index;
				// start new line for lastLeft
				fprintf(writer, "%d   ", translate[lastLeft]);
			}
			// create mcl format entries describing non-zero weights
			//   of links to lastLeft
			reader += sizeof(unsigned);
			index = *((unsigned *)reader);
			reader += sizeof(unsigned);
			float weight = *((float *)reader);
			fprintf(writer, " %d:%.3f", translate[index], weight);

			// adjance to next link for mcl input
			reader = marcher->getNextEntry();
		}
		// terminate file in standard mcl format
		fprintf(writer, " $\n)\n");
		fclose(writer);
	}
}

// use interMeans arrays to normalize ortholog weights for mcl
//   also make entries in bestInters bi-directional to
//   better support input for mcl
void OrthoTop::adjustInters(unsigned* getParent) {
	unsigned loopCount = bestInters->getEntryCount();
	//  prepare to march through all ortholog links
	BlockMarch* marcher = bestInters->getFullMarch();
	char *entry = marcher->getEntry();
	// use loopCount instead of entry == NULL
	//   since the tail bestInters will grow as more members are examined
	while (0 < loopCount) {
		// get gene coordinates for storing to mcl
		//   and genome coordinates to access proper normalization factor
		unsigned low = *((unsigned *)entry);
		unsigned lowTaxon = getParent[low];
		entry += sizeof(unsigned);
		unsigned high = *((unsigned *)entry);
		unsigned highTaxon = getParent[high];
		entry += sizeof(unsigned);
		float value = *((float *)entry);
		double adjust = interMeans[highTaxon - 1][lowTaxon];
		// normalize value for mcl
		value = (float)(value / adjust);
		*((float *)entry) = value;
		// add a link copy that goes in the opposite direction
		bestInters->add(high, low, value);
		// move to next entry
		entry = marcher->getNextEntry();
		// check to see of all original entries have been modified
		--loopCount;
	}
	// release mean arrays that are no longer required
	unsigned stop = taxonCount - 1;
	for (unsigned i = 0; i < stop; i++) {
		free(interMeans[i]);
	}
	free(interMeans);
}



// check to see if 2 genes can be linked via there paralogs
//   and a primary ortholog link
bool OrthoTop::doHubCheck(unsigned low, std::vector<unsigned> highClones){
	//  need to check primary ortholog links
	Blocks* search = bestInters->getStorage();

	// prepare minCheck to search for smallest hub pair in search
	minCheck.setToFind(low, highClones[0]);

	//  find smallest position in search that is <= pair in minCheck
	bool result = search->chopBottom(&minCheck);
	//  if minCheck was found or minCheck above all of search
	//     the work is complete
	if ((!result) && (!minCheck.aboveAll)) {
		// below will mark top of lower positions in hub where search failed
		unsigned below = 0;
		// above will mark bottom of upper position in hub where search failed
		unsigned above = highClones.size() - 1;
		// set maxCheck to search for largest hub pair
		maxCheck.setToFind(low, highClones[above]);
		//  find largest position in search that is not larger than maxCheck
		//   since minCheck tracks search position of search entries
		//      that were smaller than smallest hub entry these
		//      do not need comparison to maxCheck
		result = search->chopTop(&minCheck, &maxCheck);
		//  if minCheck found or if it is below all remaining
		//     entries in search work is complete
		if ((!result) && (!(maxCheck.belowAll))) {
			// objective now is to reduce search range bounded by
			//    minCheck & maxCheck, until match is found
			//    or range disappears
			// first try chopping search range from below
			++below;
			// loop while candidates for finding or chopping still exist
			while (below < above) {
				// set toFind to chop
				toFind.setToFind(low, highClones[below]);
				result = search->chopLower(&toFind, &minCheck, &maxCheck);
				if ((result) || (toFind.aboveAll)) {
					// if found or search range is emptied
					//   set values to exit loop
					below = above;
				} else {
					// try chopping top of search range
					--above;
					if (below <above) {
						toFind.setToFind(low, highClones[above]);
						result =
							search->chopUpper(&toFind, &minCheck, &maxCheck);
						if ((result) || (toFind.belowAll)) {
							// if found or search range is emptied
							//   set values to exit loop
							below = above;
						} else {
							// adjust below to chop lower section of
							//   range at top of loop
							++below;
						}
					}
				}
			}
		}
	}

	return result;
}

bool OrthoTop::doHubCheck(std::vector<unsigned> lowClones, unsigned high) {
	Blocks* search = bestInters->getStorage();
	// prepare minCheck to search for smallest hub pair in search
	minCheck.setToFind(lowClones[0], high);
	//  find largest position in search that is <= pair in minCheck
	bool result = search->chopBottom(&minCheck);
	//  if minCheck was found or minCheck above all of search
	//     the work is complete
	if ((!result) && (!(minCheck.aboveAll))) {
		// below will mark top of lower positions where search failed
		unsigned below = 0;
		// above will mark bottom of upper position where search failed
		unsigned above = lowClones.size() - 1;
		// set maxCheck to search for largest hub pair
		maxCheck.setToFind(lowClones[above], high);
		//  find smallest position in search that is not larger than maxCheck
		//   since minCheck tracks search position of search entries
		result = search->chopTop(&minCheck, &maxCheck);
		//  if minCheck found or if it is below all remaining
		//     entries in search work is complete
		if ((!result) && (!(maxCheck.belowAll))) {
			// objective now is to reduce search range bounded by
			//    minCheck & maxCheck, until match is found
			//    or range disappears
			// first try chopping search range from below
			++below;
			// loop while candidates for searching / chopping still exist
			while (below < above) {
				// set toFind to chop
				toFind.setToFind(lowClones[below], high);
				result = search->chopLower(&toFind, &minCheck, &maxCheck);
				if ((result) || (toFind.aboveAll)) {
					// if found or search range is emptied
					//   set values to exit loop
					below = above;
				} else {
					// try chopping top of search range
					--above;
					if (below <above) {
						toFind.setToFind(lowClones[above], high);
						result =
							search->chopUpper(&toFind, &minCheck, &maxCheck);
						if ((result) || (toFind.belowAll)) {
							// if found or search range is emptied
							//   set values to exit loop
							below = above;
						} else {
							// adjust below to chop lower section of
							//   range at top of loop
							++below;
						}
					}
				}
			}
		}
	}
	return result;
}

bool OrthoTop::doHubCheck(std::vector<unsigned> lowClones,
				std::vector<unsigned> highClones) {
	Blocks* search = bestInters->getStorage();
	// prepare minCheck to search for smallest hub pair in search
	minCheck.setToFind(lowClones[0], highClones[0]);
	bool result = search->chopBottom(&minCheck);
	//  if minCheck was found or minCheck above all of search
	//     the work is complete
	if ((!result) && (!minCheck.aboveAll)) {
		// variables for marking smallest checked member of hub
		unsigned lowBelow = 0;
		unsigned highBelow = 0;
		// support for advancing lowBelow, highBelow
		unsigned highLoop = highClones.size() - 1;

		// variables for marking largest checked member of hub
		unsigned lowAbove = lowClones.size() - 1;
		unsigned highAbove = highLoop;

		// try search for largest member of hub
		maxCheck.setToFind(lowClones[lowAbove], highClones[highAbove]);
		//  find largest position in search that is not larger than maxCheck
		result = search->chopTop(&minCheck, &maxCheck);

		//  if maxCheck found or if it is below all remaining
		//     entries in search work is complete
		if ((!result) && (!maxCheck.belowAll)) {
			// objective now is to reduce search range bounded by
			//    minCheck & maxCheck, until match is found
			//    or range disappears
			// advance to second smallest hub member
			++highBelow;
			bool done = false;
			while (!done) {
				toFind.setToFind(lowClones[lowBelow], highClones[highBelow]);
				result = search->chopLower(&toFind, &minCheck, &maxCheck);

				done = ((result) || (toFind.aboveAll));
				if (!done) {
					//  (lowBelow, highBelow) was not found and
					//      entire remaining range is not below it

					// now reduce range from above
					//    using entry below (lowAbove, highAbove)
					if (0 < highAbove) {
						--highAbove;
					} else {
						// decrement requires change most significant index
						highAbove = highLoop;
						--lowAbove;
					}
					// block check if this reaches known bottom
					//   non-match (lowBelow, highBelow)
					done = ((lowBelow == lowAbove)
								&& (highBelow == highAbove));
					if (!done) {
						// set next hub member for which to search
						toFind.setToFind(lowClones[lowAbove],
										 highClones[highAbove]);
						// search for toFind with possible range reduction
						//    if not found
						result =
							search->chopUpper(&toFind, &minCheck, &maxCheck);
						done = ((result) || (toFind.belowAll));
						if (!done) {
							// prepare (lowBelow, highBelow) for
							//    chopLower at top of loop
							if (highBelow < highLoop) {
								++highBelow;
							} else {
								// increment requires change in most
								//   significant coordinate
								highBelow = 0;
								++lowBelow;
							}
							done = ((lowBelow == lowAbove)
										&& (highBelow == highAbove));
						}
					}
				}
			}
		}
	}
	return result;
}

void OrthoTop::doOrthoConvert(const char *mclOutput, const char *idxPath,
							  const char *readable) {
	GgParser idxToName(idxPath);
	LinePull puller(0x100000, mclOutput);
	unsigned geneCount = idxToName.getGeneCount();

	// create memory for storing gene indexes found in a line
	unsigned* inLine =(unsigned *)malloc(geneCount * sizeof(unsigned));

	unsigned taxonCount = idxToName.getTaxonCount();
	// set array for tracking genomes found in a line
	char* parents = (char *)malloc(taxonCount * sizeof(unsigned));
	for (unsigned i = 0; i < taxonCount; i++) {
		parents[i] = '\0';
	}

	// skip lines header lines
	// set flag to help ignore output heading
	bool beforeBegin = true;
	char *lineStart = puller.getFrontWhiteSplit();
	while ((beforeBegin) && (lineStart != NULL)) {
		if (strcmp(lineStart, "begin") == 0) {
			beforeBegin = false;
		}
		// discard remainder of line
		while (NULL != puller.getNextWhiteSplit());
		lineStart = puller.getFrontWhiteSplit();
	}

	FILE * writer = fopen(readable, "w");
	unsigned groupCount = 0;
	while (lineStart != NULL) {
		if (*lineStart == ')') {
			// have reached end of data
			lineStart = NULL;
		} else {
			unsigned depth = 0;
			// at start of group description line
			//  group index in lineStart is not used
			char* linePart = puller.getNextWhiteSplit();
			// every line will have at least one part or no group
			while (linePart != NULL) {
				// check for end of line
				if (*linePart != '$') {
					// save index of gene in group
					inLine[depth] = atoi(linePart);
					// indicate that its parent is the lines genomes
					parents[idxToName.getParentIndex(inLine[depth])]
						        = '\1';
					++depth;
					linePart = puller.getNextWhiteSplit();
					if (linePart == NULL) {
						//  output may place genes on more than 1 line
						linePart = puller.getFrontWhiteSplit();
					}
				} else {
					// at end of ortho group
					unsigned parentCount = 0;
					for (unsigned i = 0; i < taxonCount; i++) {
						if (parents[i] != '\0') {
							++parentCount;
							parents[i] = '\0';
						}
					}
					fprintf(writer, "ORTHOMCL%d (%d genes,%d taxa):\t",
							groupCount, depth, parentCount);
					// advance for next group name
					++groupCount;
					// set names of genes in group along with
					//   their parent names
					for (unsigned i = 0; i < depth; i++) {
						fprintf(writer, " %s(%s)",
								idxToName.getGeneName(inLine[i]),
								idxToName.getParentName(inLine[i]));
					}
					fprintf(writer, "\n");
					// discard any line data after $
					linePart = puller.getNextWhiteSplit();
					while (linePart != NULL) {
						linePart = puller.getNextWhiteSplit();
					}
					lineStart = puller.getFrontWhiteSplit();
				}
			}
		}
	}
}

OrthoTop::~OrthoTop() {
}
