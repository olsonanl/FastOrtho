/*
 * GgParser.cpp
 *
 *  Created on: Apr 8, 2010
 *      Author: mscott
 */

#include "GgParser.h"
#include <algorithm>
#include <iostream>
#include "../List/Blocks.h"
#include "../List/BlockMarch.h"
#include "../OptionPass.h"
#include <ctype.h>
#include <string.h>
#include <stdlib.h>

#define PARENT_INDEX 0
#define GENE_INDEX 1
#define GENE_LENGTH 2
#define INT_COUNT 3

#define ENTRY_SIZE (sizeof(char *) + (INT_COUNT * sizeof(unsigned)))


GgParser::GgParser(const char *indexPath) {
	genomeNames = NULL;
	taxonNameRam = NULL;
	lastTaxon = NULL;
	blastWriter = NULL;

	geneLinks = NULL;
	geneNameRam = NULL;

	ggFile = NULL;
	getParent = NULL;

	indexPuller = NULL;
	indexWriter = NULL;

	geneLengths = NULL;



	// lengths are not required when loading from index file
	ggFile = indexPath;
	// associate gene names with indexes
	//   block sorting by gene name since desired gene order
	//   is order found in the index file
	loadGgEntries(false, indexPath);
}

//  constructor when loading .gg(or .glg) for use in converting
//     gene names into indices
GgParser::GgParser(const char *ggSet, const char *indexSet) {
	genomeNames = NULL;
	taxonNameRam = NULL;
	lastTaxon = NULL;
	blastWriter = NULL;

	geneLinks = NULL;
	geneNameRam = NULL;

	ggFile = NULL;
	getParent = NULL;

	indexPuller = NULL;
	indexWriter = NULL;

	geneLengths = NULL;
	ggFile = ggSet;
	indexWriter = NULL;
	indexFile = indexSet;
	// true will sort gene names to simplify conversion
	loadGgEntries(true, ggFile);
}

// preparation for absorbing gene lengths from a .bpo file
void GgParser::createLengthsRoom() {
	geneLengths = (unsigned *)malloc(geneCount * sizeof(unsigned));
}

// handles setting individual gene length
const char* GgParser::addBpoLength(const char *oldName, const char *toFind,
									const char *length) {
	// block search if length is known from last call to addBpoLength
	if ((oldName == NULL) || (0 != strcmp(oldName, toFind))) {
		GgEntry *found = getGgEntry(toFind);
		if (found != NULL) {
			geneLengths[found->order] = atoi(length);
			oldName = found->name;
		}
	}
	return oldName;
}

// handles detection of genome name in data absorbed from the end of a
//  > line from i.e. portion inside of []
char *findTaxonName(char *tail) {
	char *endChop = tail + strlen(tail) - 1;
	// discard all text which follows final ]
	while ((tail < endChop) && (*endChop != ']')) {
		--endChop;
	}
	// endChop either matches tail or points at a ]
	if (tail < endChop) {
		// discard final ]  and any _ that precede it
		--endChop;
		while ((tail < endChop) && (*endChop == '_')) {
			--endChop;
		}
	}
	if (tail == endChop) {
		*tail = '\0';
	} else {
		// set termination to create taxon name string
		*endChop = '\0';
		// have one ], looking for matching [
		//    which will be found when matchCount drops to 0
		unsigned matchCount = 1;
		while ((tail <= endChop) && (0 < matchCount)) {
			if (*endChop == ']') {
				//  puts matching [ farther to the left
				++matchCount;
			} else if (*endChop == '[') {
				//  have a match for the last ]
				--matchCount;
			}
			--endChop;
		}
		if (0 < matchCount) {
			// did not find a match for the final ]
			*tail = '0';
		} else {
			// move back to [
			++endChop;
			// move right of matching [
			++endChop;
			tail = endChop;
			// skip over whitespace at start of name
			while (*tail == '_') {
				++tail;
			}
		}
	}
	if (*tail == '\0') {
		//  send warning for empty gene names
		std::cout << "!!!!!!!" << std::endl;
		std::cout << tail << std::endl;
		// no taxon name contents
	}
	return tail;
}

// enables sorting by name when name is the first component
int compareFrontText(const char *left, const char *right) {
	const char *leftName = ((const char **)left)[0];
	const char *rightName = ((const char **)right)[0];
	return (strcmp(leftName, rightName));
}

// enables sorting genome name entries by index
int compareTaxon(const char *left, const char *right) {
	left += sizeof(char *);
	int result = ((int *)left)[0];
	right += sizeof(char *);
	result -= ((int *)right)[0];
	return result;
}

// Responsible for reading protein files
//    to get gene names, lengths, and taxon membership
// Can also be used to combine files for use as input to blast executables

//   Only 2 types of lines are recognized
//      lines that start with the > character provide gene names
//         and optionally the name of the genome which contains the gene
//     all other lines should contain the aminos in the protein sequence
void GgParser::processFasta(char *taxonNameBuild, LinePull *puller,
		Blocks *getGenes, Blocks *getTaxons) {
	char *next = puller->getFrontWhiteSplit();
	char *geneName = NULL;
	unsigned geneLength = 0;
	while (next != NULL) {
		if (*next != '>') {
			// use non > lines to find gene length
			do {
				geneLength += strlen(next);
				// put protein string in file for blast executables
				if ((geneName != NULL) && (blastWriter != NULL)) {
					fprintf(blastWriter, "%s\n", next);
				}
				next = puller->getNextWhiteSplit();
			} while (next != NULL);
		} else {
			// at gene & taxon name line
			//  check to see if data needs to be stored for previous gene
			if ((geneName != NULL) && (0 < geneLength)) {
				// add previous gene data to getGenes
				char* storeGene = getGenes->getAdd();
				//  name of gene
				*((const char **)storeGene) = geneName;
				storeGene += sizeof(char *);

				// store information about gene
				unsigned *intStore = ((unsigned *)storeGene);
				// store taxon membership of gene
				intStore[PARENT_INDEX] = taxonCount - 1;
				//  store index of gene
				intStore[GENE_INDEX] = geneCount;
				//  store sequence length of gene
				intStore[GENE_LENGTH] = geneLength;
				// increase gene count for this gene
				++geneCount;
			}
			// check on placing line start in combined blast file
			if (blastWriter != NULL) {
				fprintf(blastWriter, "%s", next);
			}
			// gene name is terminated at whitespace but does not
			//   	include > intro
			// gene name needs to be in ram that will not be modified
			//     unlike the memory at next
			geneName = geneNameRam->add(next + 1);
			// reset gene length
			geneLength = 0;

			// deal with remainder of line
			next = puller->getNextWhiteSplit();
			if (taxonNameBuild == NULL) {
				// do not need to extract taxon name
				//  so just check on passing parts to combined blast file
				//     multi-character whitespace will be replaced by
				//     a single space character
				while (next != NULL) {
					if (blastWriter != NULL) {
						fprintf(blastWriter, " %s", next);
					}
					next = puller->getNextWhiteSplit();
				}
			} else {
				//  use buffer to construct taxon name
				char *buildAdd = taxonNameBuild;
				// place remainder of line in a single string
				if (next == NULL) {
					// if nothing in line taxon name will be empty
					printf("No genome for %s\n", geneName);
					buildAdd[0] = '\0';
				} else {
					// join remainder of line into a single string
					//   replacing white space blocks with '_'
					do {
						// pass remainder of line to combined blast file
						if (blastWriter != NULL) {
							fprintf(blastWriter, " %s", next);
						}
						// enlarge string from which
						//    taxon name will be extracted
						unsigned span = strlen(next);
						memcpy(buildAdd, next, span);
						buildAdd += span;
						// add to separate from next non-white block
						*buildAdd = '_';
						++buildAdd;
						next = puller->getNextWhiteSplit();
					} while (next != NULL);
					// drop trailing _
					--buildAdd;
					*buildAdd = '\0';
					//  find text that is enclosed by teminating []
					buildAdd = findTaxonName(taxonNameBuild);
				}
				// check for change in taxon name
				if ((lastTaxon == NULL) ||
						(0 != strcmp(buildAdd, lastTaxon))) {
					// get a copy of taxon name in ram
					// that will not be modified
					lastTaxon = taxonNameRam->add(buildAdd);

					// save new taxon name and its discovery position
					char *storeTaxon = getTaxons->getAdd();
					*((const char **)storeTaxon) = lastTaxon;
					storeTaxon += sizeof(char *);

					*((unsigned *)storeTaxon) = taxonCount;
					++taxonCount;
				}
			}
			// terminate gene name line in combined fasta file
			if (blastWriter != NULL) {
				fprintf(blastWriter, "\n");
			}
		}
		// move on to next protein file line
		next = puller->getFrontWhiteSplit();
	}
	// check on saving data about final gene
	if ((geneName != NULL) && (0 < geneLength)) {
		// add previous gene data to getGenes
		char* storeGene = getGenes->getAdd();
		//  name of gene
		*((char **)storeGene) = geneName;

		storeGene += sizeof(char *);
		unsigned *intStore = ((unsigned *)storeGene);
		intStore[PARENT_INDEX] = taxonCount - 1;
		intStore[GENE_INDEX] = geneCount;
		intStore[GENE_LENGTH] = geneLength;
		++geneCount;
	}
}

// construct that starts with a collection of proteing files
GgParser::GgParser(unsigned singleCount, const char **genomeFastas,
					 unsigned mixedCount, const char **mixedFastas,
					 const char *ggStore, const char *indexSet,
					 const char *combined) {
	genomeNames = NULL;
	lastTaxon = NULL;

	geneLinks = NULL;

	getParent = NULL;

	indexPuller = NULL;
	indexWriter = NULL;

	geneLengths = NULL;



	ggFile = ggStore;
	indexFile = indexSet;
	// check to see if combining all fasta files into a single
	//   file has been requested
	if (combined == NULL) {
		blastWriter = NULL;
	} else {
		blastWriter = fopen(combined, "w");
	}

	// create ram for capture text of taxon names and gene names
	taxonNameRam = new StringStore(0x10000);
	geneNameRam = new StringStore(0x10000);

	// in getGenes store name, parent index, gene index, gene length
	Blocks getGenes(ENTRY_SIZE);

	// for each taxon store order of discovery (parent index in gene)
	//    to handle possibility that all genes for a single
	//    taxon are not stored in a single block
	//
	Blocks getTaxons(sizeof(char*) + sizeof(unsigned));

	taxonCount = 0;
	lastTaxon = NULL;

	// first handle case where taxon names are derived from
	//    file names
	geneCount = 0;
	for (unsigned i = 0; i < singleCount; i++) {
		const char *fastaFile = genomeFastas[i];
		const char *nameStop = fastaFile + strlen(fastaFile);
		const char *start = nameStop - 1;
		// remove directory path from name
		while ((fastaFile < start) && (*start != '/') && (*start != '\\')) {
			--start;
		}
		if (fastaFile < start) {
			// skip over directory boundary
			++start;
		}
		// remove file tail from name i.e. drop ".faa" at end
		const char *check = nameStop - 1;
		while (start < check) {
			if (*check != '.') {
				--check;
			} else {
				nameStop = check;
				check = start;
			}
		}
		// put constructed name in ram that will not get modified
		lastTaxon =  taxonNameRam->add(start, nameStop);
		//  add taxon information to getTaxons
		char *storeTaxon = getTaxons.getAdd();
		*((char **)storeTaxon) = lastTaxon;
		storeTaxon += sizeof(char *);

		*((unsigned *)storeTaxon) = taxonCount;
		++taxonCount;
		//  extract gene names and sequences from fasta file
		fastaFile = findFasta(fastaFile);
		LinePull puller(0x10000, fastaFile);
		processFasta(NULL, &puller, &getGenes, &getTaxons);
	}

	// create work buffer for isolating taxon name
	//   from contents of a > line
	char* taxonNameBuild = (char *)malloc(0x1000 * sizeof(char));
	for (unsigned i = 0; i < mixedCount; i++) {
		const char *fastaFile = findFasta(mixedFastas[i]);
		LinePull puller(0x10000, fastaFile);
		processFasta(taxonNameBuild, &puller, &getGenes, &getTaxons);
	}
	// close combined file if it was being created
	if (blastWriter != NULL) {
		fclose(blastWriter);
	}
	//  release buffer used to build taxon names from > lines
	free(taxonNameBuild);
	// remove excess from ram that stores taxon & genome names
	taxonNameRam->releaseTail();
	geneNameRam->releaseTail();

	// all of the taxon names are stored in getTaxons, but it is possible
	//   for a single name to appear multiple times in getTaxons
	// To check for this difficulty the code will sort the getTaxon
	//    entries by the taxon names and then make adjustments
	//    to insure that genes linked to the same taxon name will
	//    also be linked to the same taxon index

	//  In getTaxons there is a value that links each taxon name
	//    to a block index of a group of identical taxon names.
	//  The taxonMap will provide a means of converting these block
	//     appearance indices to unique taxon alphabetic indices
	unsigned *taxonMap = (unsigned *)malloc(taxonCount * sizeof(unsigned));
	genomeNames = (char **)malloc(taxonCount * sizeof(char *));

	// sort getTaxons alphabetically based on  taxon names
	getTaxons.sort(compareFrontText);

	// march through getTaxon to set taxonMap and to fill genomeNames
	BlockMarch* marcher = getTaxons.getFullMarch();
	char *reader = marcher->getEntry();
	// get first taxon name (based on alphabetic order)
	char *lastName = *((char **)reader);
	genomeNames[0] = lastName;
	//  set first entry in taxonMap
	reader += sizeof(char *);
	taxonMap[*((unsigned *)reader)] = 0;

	//  track number of non-identical taxon names
	unsigned uniqueTaxon = 0;
	reader = marcher->getNextEntry();
	while (reader != NULL) {
		char *nextName = *((char **)reader);
		// check for a change in taxon name
		if (strcmp(lastName, nextName) != 0) {
			// update last name
			lastName = nextName;
			// add new entry to genomeNames
			++uniqueTaxon;
			genomeNames[uniqueTaxon] = lastName;
		}
		reader += sizeof(char *);
		//  map appearance of a taxon block to alphabetic order
		//     of name for block
		taxonMap[*((unsigned *)reader)] = uniqueTaxon;
		reader = marcher->getNextEntry();
	}
	// convert from last storage position to entries in genomeNames
	++uniqueTaxon;

	// The stored gene information includes values that link
	//    the gene to taxon name blocks which shall be replaced
	//    by a value that describes the alphabetic position of its associated
	//    gene name.
	marcher = getGenes.getFullMarch();
	reader = marcher->getEntry();
	while (reader != NULL) {
		// skip over pointer to gene name
		reader += sizeof(char *);
		// get taxon name block location
		unsigned index = *((unsigned *)reader);
		// convert to alphabetic order for taxon name
		index = taxonMap[index];
		//  save conversion
		*((unsigned *)reader) = index;
		reader = marcher->getNextEntry();
	}
	//  discard map
	free(taxonMap);

	if (uniqueTaxon < taxonCount) {
		// protein files produced non-adjacent blocks of genome names
		//  with identical taxon names so reduce taxonCount
		//    and release unused memory in genomeNames
		taxonCount = uniqueTaxon;
		genomeNames =
				(char **)realloc(genomeNames, taxonCount * sizeof(char *));
	}
	// sort genes according to taxon index to insure
	//     all genes for a single taxon appear in a single block
	getGenes.sort(compareTaxon);

	// replace discovery position stored in getGene
	//    with a position that is based on its genome
	marcher = getGenes.getFullMarch();
	reader = marcher->getEntry();
	unsigned newOrder = 0;
	while (reader != NULL) {
		// skip over name pointer and parent order
		reader += sizeof(char *) + sizeof(unsigned);
		*((unsigned *)reader) = newOrder;
		++newOrder;
		reader = marcher->getNextEntry();
	}
	// create a file that will allow reloading of gene names
	//    and their associated genome names and sequence lengths.
	createGgFile(&getGenes);
}

const char *testName(const char *toCheck) {
	FILE* tester = fopen(toCheck, "r");
	if (tester == NULL) {
		toCheck = NULL;
	} else {
		fclose(tester);
	}
	return toCheck;
}

const char *GgParser::findFasta(const char *name) {
	const char *result = testName(name);
	if (result == NULL) {
		unsigned baseSpan = strlen(name);
		// create enough room to add .fasta to end of name
		char *buffer = (char *)malloc(7 + baseSpan);
		memcpy(buffer, name, baseSpan);
		strcpy(buffer + baseSpan, ".fasta");
		result = testName(buffer);
		if (result == NULL) {
			strcpy(buffer + baseSpan, ".fa");
			result = testName(buffer);
			if (result == NULL) {
				strcpy(buffer + baseSpan, ".faa");
				result = testName(buffer);
			} else {
				std::cerr << "Unable to find " << name << std::endl;
				strcpy(buffer + baseSpan, ".fasta");
				std::cerr << " or " << buffer << std::endl;
				strcpy(buffer + baseSpan, ".fa");
				std::cerr << " or " << buffer << std::endl;
				strcpy(buffer + baseSpan, ".faa");
				std::cerr << " or " << buffer << std::endl;
				free(buffer);
				// This function generates a memory leak if
				//    buffer is used for the result
				//    but the amount of memory lost should not
				//    be significant
			}
		}
	}
	return result;
}

void GgParser::createGgFile(Blocks *getGenes) {
	// create a file to save this information for reload
	FILE * writer = fopen(ggFile, "w");
	// set lastTaxon so the first discovered taxon
	//    will be seen as new by lastTaxon will be seen as invalid
	unsigned lastTaxon = taxonCount;
	unsigned geneAt = 0;
	// create array for quick access to gene lengths
	geneLengths = (unsigned *)malloc(geneCount * sizeof(unsigned));

	BlockMarch * marcher = getGenes->getFullMarch();
	char *reader = marcher->getEntry();
	while (reader != NULL) {
		char *geneName = *((char **)reader);
		reader += sizeof(char *);
		unsigned *intRead = (unsigned *)reader;
		unsigned nextTaxon = intRead[PARENT_INDEX];
		if (nextTaxon != lastTaxon) {
			// need to end previous line if valid
			if (lastTaxon < taxonCount) {
				fprintf(writer, "\n");
			}
			//  start new taxon line
			lastTaxon = nextTaxon;
			fprintf(writer, "%s:", genomeNames[lastTaxon]);
		}
		// add gene length and name
		fprintf(writer, " {%d}%s", intRead[GENE_LENGTH], geneName);
		geneLengths[geneAt] = intRead[GENE_LENGTH];

		++geneAt;
		reader = marcher->getNextEntry();
	}
	fprintf(writer, "\n");
	fclose(writer);

	//  sort getGenes to aid search for gene information based on
	//     gene name
	getGenes->sort(compareFrontText);

	// create GgEntry array that will be used to get
	//    gene information from a gene name
	marcher = getGenes->getFullMarch();
	reader = marcher->getEntry();
	geneLinks = (GgEntry *)(malloc(geneCount * sizeof(GgEntry)));
	for (unsigned i = 0; i < geneCount; i++) {
		//  first copy gene name
		(geneLinks[i]).name = *((char **)reader);
		// next copy genome and gene indices for name
		unsigned* intRead = (unsigned *)(reader + sizeof(char *));
		(geneLinks[i]).parent = intRead[PARENT_INDEX];
		(geneLinks[i]).order = intRead[GENE_INDEX];
		reader = marcher->getNextEntry();
	}
}

void GgParser::getGgIndexes() {
	loadGgEntries(false, indexFile);
}

unsigned GgParser::getParentIndex(unsigned geneIndex) {
	return ((geneLinks[geneIndex]).parent);
}

const char * GgParser::getParentName(unsigned geneIndex) {
	unsigned index = (geneLinks[geneIndex]).parent;
	return (genomeNames[index]);
}
const char * GgParser::getGeneName(unsigned geneIndex) {
	return ((geneLinks[geneIndex]).name);
}

// creates geneLinks from special files that describe
//    relationship between genes and genomes
void GgParser::loadGgEntries(bool doSort, const char *ggPath) {
	bool haveGeneLengths = true;
	// create areas for storing taxon name and
	//    geneome name text discovered in gg file
	taxonNameRam = new StringStore(0x10000);
	geneNameRam = new StringStore(0x10000);

	Blocks taxonList(sizeof(char *));

	//  create storage for information about genome names
	Blocks collect(ENTRY_SIZE);

	//  get access to .gg  or .glg file
	LinePull puller(0x10000, ggPath);

	taxonCount = 0;
	geneCount = 0;
	char *next = puller.getFrontWhiteSplit();
	while (next != NULL) {
		// first entry in line might be a genome name
		int tail = strlen(next) - 1;
		if (next[tail] == ':') {
			// first entry in line is a genome name
			//  discard :  taxon name indicator
			next[tail] = '\0';
			char *addAt = taxonList.getAdd();
			// save address of unchangable copy of taxon name from file
			((char **)addAt)[0] = taxonNameRam->add(next);
			// update number of taxon names detectect
			++taxonCount;
			//  set next for finding genome names
			next = puller.getNextWhiteSplit();
		}
		// remainder of line entries must be genes
		while (next != NULL) {
			// do no save gene information if
			//   it is does not follow a taxon name
			if (0 < taxonCount) {
					unsigned geneLength = 0;
				if (*next != '{') {
					haveGeneLengths = false;
				} else {
					char *digiter = next + 1;
					while (isdigit(*digiter)) {
						++digiter;
					}
					if ((*digiter != '}') || (next + 1 == digiter)) {
						// either digits did not end with }
						//  or prefix looked like {}
						haveGeneLengths = false;
					} else {
						*digiter = '\0';
						// get gene length from {length} prefix to gene name
						geneLength = atoi(next + 1);
						// set next to ignore prefix
						next = digiter + 1;
					}
				}
				// get block location for storing gene information
				char *store = collect.getAdd();
				// store name address in geneNameRoom
				char *nameAt = NULL;
				nameAt = geneNameRam->add(next);

				*((char **)store) = nameAt;
				store += sizeof(char *);

				unsigned * intStore = (unsigned *)store;
				intStore[PARENT_INDEX] = taxonCount - 1;

				intStore[GENE_INDEX] = geneCount;
				++geneCount;

				intStore[GENE_LENGTH] = geneLength;
			}
			next = puller.getNextWhiteSplit();
		}
		next = puller.getFrontWhiteSplit();
	}
	geneNameRam->releaseTail();
	taxonNameRam->releaseTail();
	printf("gene count = %d in %d taxons", geneCount, taxonCount);
	std::cout << std::endl;

	genomeNames = (char **)(malloc(taxonCount * sizeof(char *)));
	BlockMarch *marcher = taxonList.getFullMarch();
	char *reader = marcher->getEntry();
	for(unsigned i = 0; i < taxonCount; i++) {
		genomeNames[i] = ((char **)reader)[0];
		reader = marcher->getNextEntry();
	}


	// If entries are being collected to convert output of mcl
	//   processing desired order for geneLinks is order provided
	//   in file. (doSort == false)
	// However, if entries are to support conversion of gene names
	//   into indices prior to mcl process order of geneLinks
	//   should be base on alphabetic order of gene names (doSort == true)
	if (doSort) {
		collect.sort(compareFrontText);
	}

	// collect is now ready to create geneLinks
	marcher = collect.getFullMarch();
	reader = marcher->getEntry();
	geneLinks = (GgEntry *)(malloc(geneCount * sizeof(GgEntry)));
	if (haveGeneLengths) {
		geneLengths = (unsigned *)(malloc(geneCount * sizeof(unsigned)));
		for (unsigned storeAt = 0; storeAt < geneCount; storeAt++) {
			(geneLinks[storeAt]).name = *((char **)reader);
			reader += sizeof(char *);

			unsigned *intRead = (unsigned *)reader;

			(geneLinks[storeAt]).parent = intRead[PARENT_INDEX];
			(geneLinks[storeAt]).order = intRead[GENE_INDEX];
			geneLengths[intRead[GENE_INDEX]] = intRead[GENE_LENGTH];

			reader = marcher->getNextEntry();
		}

	} else {
		for (unsigned storeAt = 0; storeAt < geneCount; storeAt++) {
			(geneLinks[storeAt]).name = *((char **)reader);
			reader += sizeof(char *);

			unsigned *intRead = (unsigned *)reader;

			(geneLinks[storeAt]).parent = intRead[PARENT_INDEX];
			(geneLinks[storeAt]).order = intRead[GENE_INDEX];
			reader = marcher->getNextEntry();
		}
	}
}


void GgParser::addToIndexFile(unsigned geneIndex, unsigned mclIndex) {
	if (indexWriter == NULL) {
		// adding first gene to index file
		//   it will be necessary to use original .gg or .glg file
		//   to extract names to place in the index file
		indexPuller = new LinePull(0x100000, ggFile);
		taxonName = strdup(indexPuller->getFrontWhiteSplit());
		indexAt = 0;
		// create file for saving .idx for conversion after mcl processing
		indexWriter = fopen(indexFile, "w");
		haveGenesForTaxon = false;
	}
	// skip through indexPuller entries until we reach geneIndex
	while (indexAt < geneIndex) {
		char *nextGeneName = indexPuller->getNextWhiteSplit();
		if (nextGeneName == NULL) {
			// found all entries for taxonName
			//    check to see if the end of its line needs to be set
			if (haveGenesForTaxon) {
				fprintf(indexWriter, "\n");
				haveGenesForTaxon = false;
			}
			free(taxonName);
			// set next taxonName
			taxonName = strdup(indexPuller->getFrontWhiteSplit());
		} else {
			// skip over nextGeneName
			++indexAt;
		}
	}
	// next gene from indexPuller goes with geneIndex
	//   but we might still need to flush an old taxon
	char *writeName = indexPuller->getNextWhiteSplit();
	while (writeName == NULL) {
		//  set about flushing output for taxonName
		if (haveGenesForTaxon) {
			fprintf(indexWriter, "\n");
			haveGenesForTaxon = false;
		}
		free(taxonName);
		// load next taxonName
		taxonName = strdup(indexPuller->getFrontWhiteSplit());
		// and its first gene
		writeName = indexPuller->getNextWhiteSplit();
	}
	// check to see if this is the first gene for taxonName
	if (!haveGenesForTaxon) {
		// note: taxonName gets it : ending from indexPuller
		fprintf(indexWriter, "%s", taxonName);
		haveGenesForTaxon = true;
	}
	// append gene name to end of line
	fprintf(indexWriter, " %s", writeName);
	++indexAt;
}

void GgParser::closeIndexFile() {
	if (taxonName != NULL) {
		free(taxonName);
	}
	if (haveGenesForTaxon) {
		fprintf(indexWriter, "\n");
	}
	fclose(indexWriter);
	indexWriter = NULL;
}


// provides information about a gene name
GgEntry* GgParser::getGgEntry(const char * geneName) {
	//  set result with assumption that name will not be found
	GgEntry* result = NULL;
	if (0 < geneCount) {
		unsigned above = geneCount - 1;
		//  first check to see if geneName matches
		//   largest geneLinks name (last alphabetically)
		int differ = strcmp(geneName, (geneLinks[above]).name);
		// 0 < differ means above all geneLinks so NULL result is correct
		if (differ <= 0) {
			if (differ == 0) {
				// return entry matching geneName
				result = &(geneLinks[above]);
			} else {
				// next check to see if geneName matches smallest
				//  geneLinks name (1st alphabetically)
				differ = strcmp(geneName, (geneLinks[0]).name);
				// differ < 0 means below all geneLinks so keep NULL result
				if (0 <= differ) {
					if (differ == 0) {
						// return entry matching geneName
						result = &(geneLinks[0]);
					} else {
						int below = 0;
						// have geneLinks entries below and above
						//   geneName
						// These bounds will be reduced in a binary
						//   fashion until and entry match for geneName
						//   is found it is proven that such a match
						//    can not happen
						// Split interval that might contain match
						int check = (below + above) >> 1;
						//  work is over when interval will not split
						//    or if match is discovered
						while (check != below) {
							// determine which side of the middle any
							//   match must be found
							differ = strcmp(geneName, (geneLinks[check]).name);
							if (differ < 0) {
								// match can only be in lower half of interval
								above = check;
								// set middle for lower half of interval
								check = (below + above) >> 1;
							} else if (differ == 0) {
								// found match
								result = &(geneLinks[check]);
								//  set to exit loop
								check = below;
							} else {
								// match can only be in upper half of interval
								below = check;
								// set middle for upper half of interval
								check = (below + above) >> 1;
							}
						}
					}
				}
			}
		}
	}
	if (result == NULL) {
		std::cerr << geneName << " is not defined in .gg or .glg file" << std::endl;
	}
	return result;
}


// after m8 or bpo file has been parsed
//    the ram which deals with the text version of taxon names or
//    genome names is no longer required and can be released
void GgParser::releaseNameRam() {
	free(genomeNames);
	genomeNames = NULL;
	delete(taxonNameRam);
	taxonNameRam = NULL;

	// save index link from gene index to taxon index
	getParent = (unsigned *)(malloc(geneCount * sizeof(unsigned)));
	for (unsigned i = 0; i < geneCount; i++) {
		GgEntry* next = &(geneLinks[i]);
		getParent[next->order] = next->parent;
	}
	free(geneLinks);
	geneLinks = NULL;

	delete(geneNameRam);
	geneNameRam = NULL;
}

void GgParser::releaseParentIndex() {
	free(getParent);
	getParent = NULL;
}

unsigned GgParser::getTaxonCount() {
	return taxonCount;
}

unsigned GgParser::getGeneCount() {
	return geneCount;
}

unsigned* GgParser::getParentIndex() {
	return getParent;
}

GgParser::~GgParser() {
	if (genomeNames != NULL) {
		free(genomeNames);
		genomeNames = NULL;
	}
	if (taxonNameRam != NULL) {
		delete(taxonNameRam);
		taxonNameRam = NULL;
	}
	if (geneLinks != NULL) {
		free(geneLinks);
		geneLinks = NULL;
	}
	if (geneNameRam != NULL) {
		delete(geneNameRam);
		geneNameRam = NULL;
	}

	if (getParent != NULL) {
		free(getParent);
		getParent = NULL;
	}
}
