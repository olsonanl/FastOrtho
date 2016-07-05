/*
 * GgParser.h
 *
 *  Created on: Apr 8, 2010
 *      Author: mscott
 */

#ifndef GGPARSER_H_
#define GGPARSER_H_
#include <string>
#include "../List/GgEntry.h"
#include "../FileIO/LinePull.h"
#include "../List/StringStore.h"
#include "../List/Blocks.h"

class GgParser {
private:
	char **genomeNames;
	StringStore *taxonNameRam;
	unsigned taxonCount;
	char *lastTaxon;
	FILE *blastWriter;

	GgEntry *geneLinks;
	unsigned geneCount;

	StringStore *geneNameRam;

	const char *ggFile;
	unsigned* getParent;

	LinePull* indexPuller;
	FILE *indexWriter;
	unsigned indexAt;
	char *taxonName;
	bool haveGenesForTaxon;
	void loadGgEntries(bool doSort, const char *ggFile);
	void loadFromFasta();
	void createGgFile(Blocks *getGenes);

	void processFasta(char *taxonNameBuild, LinePull *puller,
					Blocks *getGenes, Blocks *getTaxons);

public:
	GgParser(const char *indexPath);
	GgParser(const char *ggSet, const char *indexSet);

	GgParser(unsigned singleCount, const char **genomeFastas,
			 unsigned mixedCount, const char **mixedFastas,
			 const char *ggStore, const char *indexSet,
			 const char *combined);


	void createLengthsRoom();
	const char *addBpoLength(const char *lastName, const char *geneName,
							 const char *length);

	const char *findFasta(const char *name);

	const char *indexFile;
	unsigned *geneLengths;
	GgEntry* getGgEntry(const char *geneName);
	void releaseNameRam();
	void releaseParentIndex();
	unsigned* getParentIndex();
	unsigned getParentIndex(unsigned geneIndex);
	const char * getParentName(unsigned geneIndex);
	const char * getGeneName(unsigned geneIndex);
	unsigned getTaxonCount();
	unsigned getGeneCount();
	void getGgIndexes();
	void addToIndexFile(unsigned geneIndex, unsigned mclIndex);
	void closeIndexFile();
	virtual ~GgParser();
};

#endif /* GGPARSER_H_ */
