/*
 * OptionPass.h
 *
 *  Created on: Apr 8, 2010
 *      Author: mscott
 */

#ifndef OPTIONPASS_H_
#define OPTIONPASS_H_

#include "LinkParser/LinkParser.h"
#include "Parser/GgParser.h"
#include "Executer.h"
#include "List/Links.h"
#include <stdio.h>
#define PATH_BREAK '/'

#if WIN32
#define MKDIR mkdir(result)
#else
#define MKDIR mkdir(result, 0755)
#endif

class OptionPass {
private:
	const char* optionFile;
	bool gotOrthoSetting;
	const char* pCutOpt;
	const char* piCutOpt;
	const char* matchCutOpt;
	const char* maxWeightOpt;
	const char* linkOpt;
	const char* bpoOpt;
	const char* ggOpt;
	const char* workingOpt;
	const char* projectOpt;
	const char* formatDbPath;
	const char* blastAllPath;
	const char* mclPath;
	bool setResultPath;
	const char* resultOpt;
	const char* queryIndex;
	const char* subjectIndex;
	const char* eValueIndex;
	const char* identityIndex;
	const char* areaIndex;
	const char* queryStartIndex;
	const char* queryEndIndex;
	const char* queryLengthIndex;
	const char* subjectStartIndex;
	const char* subjectEndIndex;
	const char* subjectLengthIndex;
	const char* mapColumnIndex;
	const char* tabSplit;
	const char* orthoQuirkOpt;
	const char* addGenome;
	const char* setAll;
	const char* inflationOpt;
	bool inflationSet;
	const char* cpuOpt;
	const char* bOpt;
	const char* eOpt;
	const char* vOpt;
	const char* splitOpt;
	const char* preBlastOpt;
	const char* legacyOpt;

	unsigned optFileDepth;

	void getOptionsFromFile(const char *filePath);
	bool checkOptionFile(const char *filePath);
	unsigned useOption(const char* optionType, const char *optionValue);
	char *projectFileStart;
	int projectFileFront;
	char *getProjectFile(const char *suffix);
	Executer *executer;

	void doQuirkSplit(LinkParser *parser, Links** savers);
	bool createGgFromFastas(bool result);
	bool startBlastWorkflow(bool result);
	unsigned *indexSettings;

	FILE* getOptionFile();
	FILE* storeFastaLists();
	void saveLinkDetails(FILE* optSave);
	void saveSettings(FILE* optSave);

public:
	GgParser *ggParser;
	const char * eText;
	double maxLogE;
	double maxWeight;

	unsigned mapColumn;

	double matchCut;

	double identityCut;

	unsigned singleCount;
	const char **genomeFastas;
	unsigned mixedCount;
	const char **mixedFastas;
	const char *allFasta;

	const char *ggFile;
	const char* linkFile;
	const char* workingDir;
	const char* projectName;
	const char* resultFile;

	const char* formatdb;
	const char* blastall;
	const char* mcl;

	char splitter;

	bool useQuirks;
	bool bpoInput;

	const char *inflation;

	const char *cpuCount;
	const char *bBlast;
	const char *vBlast;
	const char *eBlast;

	bool blastPrep;
	bool useLegacy;

	GgParser *getGgParser();
	char *idxFile;

	bool classifyLinks(Links **bins);

	char *getMatrixPath();
	char *matrixPath;
	char *getMclOutputPath();
	char *mclOutputPath;

	bool runMcl();

	OptionPass(int argc, const char **argv);
	bool areOk();
	virtual ~OptionPass();
};

#endif /* OPTIONPASS_H_ */
