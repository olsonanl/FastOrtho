/*
 * OptionPass.cpp
 *
 *  Created on: Apr 8, 2010
 *      Author: mscott
 */

#include "OptionPass.h"
#include <iostream>
#include <sstream>
#include <string>
#include <cstdlib>
#include <math.h>
#include "Parser/CharSplitParser.h"
#include "Parser/BpoAllParser.h"
#include "Classify/Ortho.h"
#include "Classify/Normal.h"
#include "LinkParser/PickE.h"
#include "LinkParser/PickPiE.h"
#include "LinkParser/PickAll.h"

#if WIN32
#include <dir.h>
#define STAT64 stat
#else
#define STAT64 stat64
#endif
#include <dirent.h>

#define _FILE_OFFSET_BITS 64
#include <sys/stat.h>
#include <float.h>
#include <string.h>
#include "Parser/LineParser.h"
#include <time.h>



using namespace std;

// supports get log 10 of text in possible exponential format
float getLog10(const char *eText) {
	float result = -FLT_MAX;

	const char *eSplit = strchr(eText, 'e');
	if (eSplit == NULL) {
		eSplit = strchr(eText, 'E');
	}
	if (eSplit == NULL) {
		// no need to split
		double aid = atof(eText);
		// use default return if aid has no log
		if (0.0 < aid) {
			aid = log10(aid);
			result = (float)aid;
		}
	} else {
		// put front part of eText into a string
		++eSplit;
		int preSize = eSplit - eText;
		// create room for start include e
		char *preAid = (char *)malloc(preSize);
		// copy front excluding e
		memcpy(preAid, eText, preSize - 1);
		// terminate where e would be
		preAid[preSize] = '\0';
		// extract value from string
		double aid = atof(preAid);
		free(preAid);
		// use default if value in aid does not have a log
		if (0.0 < aid) {
			aid = log10(aid);
			aid += atof(eSplit);
			result = (float)aid;
		}
	}
	return result;
}


OptionPass::OptionPass(int argc, const char **argv) {
	// initialize option descriptions
	optionFile = "--option_file";
	gotOrthoSetting = false;
	pCutOpt = "--pv_cutoff";
	piCutOpt = "--pi_cutoff";
	matchCutOpt = "--pmatch_cutoff";
	maxWeightOpt = "--maximum_weight";
	linkOpt = "--blast_file";
	bpoOpt = "--bpo_file";
	ggOpt = "--gg_file";
	workingOpt = "--working_directory";
	projectOpt = "--project_name";
	formatDbPath = "--formatdb_path";
	blastAllPath = "--blastall_path";
	mclPath = "--mcl_path";
	setResultPath = false;
	resultOpt = "--result_file";
	queryIndex = "--query_index";
	subjectIndex = "--subject_index";
	eValueIndex = "--e_value_index";
	identityIndex = "--percent_idenity_index";
	areaIndex = "--alignment_length_index";
	queryStartIndex = "--query_start_index";
	queryEndIndex = "--query_end_index";
	queryLengthIndex = "--query_length_index";
	subjectStartIndex = "--subject_start_index";
	subjectEndIndex = "--subject_end_index";
	subjectLengthIndex = "--subject_length_index";
	mapColumnIndex = "--mapping_index";
	tabSplit = "--use_tab_split";
	orthoQuirkOpt = "--match_OrthoMcl";
	addGenome = "--single_genome_fasta";
	setAll = "--mixed_genome_fasta";
	inflationSet = false;
	inflationOpt = "--inflation";
	cpuOpt = "--blast_cpus";
	bOpt = "--blast_b";
	eOpt = "--blast_e";
	vOpt = "--blast_v";
	splitOpt = "--split_char";
	preBlastOpt = "--only_fastas";
	legacyOpt = "--legacy_blast";

	// assume preBlastOpt will not appear
	blastPrep = false;
	useLegacy = false;

	// create object for making system calls
	executer = new Executer();
	matrixPath = NULL;
	mclOutputPath = NULL;

	// set blast defaults
	cpuCount = "1";
	bBlast = "1000";
	vBlast = "1000";
	eBlast = "1e-5";

	// set indexSettings with defaults for -m 8 file
	//  (and those specific for a bpo file with bpo defaults
	indexSettings = (unsigned *)malloc(INDEX_ROOM * sizeof(unsigned));
	indexSettings[QUERY_AT] = 0;
	indexSettings[QUERY_START] = 6;
	indexSettings[QUERY_END] = 7;
	indexSettings[QUERY_LENGTH] = 2;
	indexSettings[SUBJECT_AT] = 1;
	indexSettings[SUBJECT_START] = 8;
	indexSettings[SUBJECT_END] = 9;
	indexSettings[SUBJECT_LENGTH] = 4;
	indexSettings[E_VALUE_AT] = 10;
	indexSettings[IDENTITY_AT] = 2;
	indexSettings[LENGTH_AT] = 3;
	indexSettings[MAP_COLUMN] = 7;

	// set default column spliter for -m 8 file
	splitter = '\t';

	// set default maximum log of e-value
	eText = "1e-5";
	maxLogE = getLog10(eText);
	// set default value to use on e-values that do not have a valid log
	maxWeight = 316.0;

	// default is not to filter on percent identity or percen match
	identityCut = 0.0;
	matchCut = 0.0;

	// there are no default protein files
	singleCount = 0;
	genomeFastas = NULL;
	mixedCount = 0;
	mixedFastas = NULL;

	//  there are no default files
	ggFile = NULL;
	linkFile = NULL;
	workingDir = NULL;
	projectName = NULL;
	resultFile = NULL;

	// set default starting point for blast executables or mcl
	formatdb = NULL;
	blastall = NULL;
	mcl = "mcl";
	// set only option value for mcl
	inflation = "1.5";

	//  assume flow will use a -m 8 file
	bpoInput = false;

	// assume orthmcl odities will not be used
	useQuirks = false;

	optFileDepth = 0;
	if (2 < argc) {
		// march through command line options
		const char **lastArg = argv + argc;
		const char **fromLine = argv + 1;
		while (fromLine < lastArg) {
			const char *nextOpt = *fromLine;
			++fromLine;
			if (0 < strlen(nextOpt)) {
				if (fromLine < lastArg) {
					fromLine += useOption(nextOpt, *fromLine);
				} else {
					// no value for nextOpt but nextOpt may not require a value
					useOption(nextOpt, "");
				}
			}
		}
	} else if (2 == argc) {
		checkOptionFile(argv[1]);
	}
	if (formatdb == NULL) {
		if (useLegacy) {
			formatdb = "formatdb";
		} else {
			formatdb = "makeblastdb";
		}
	}
	if (blastall == NULL) {
		if (useLegacy) {
			blastall = "blastall";
		} else {
			blastall = "blastp";
		}

	}
}

bool OptionPass::checkOptionFile(const char *filePath) {
	LinePull puller(0x10000, filePath);
	bool result = false;
	unsigned lineCount = 0;
	char* nextLine = puller.getLine();
	while (nextLine != NULL) {
		nextLine = puller.trim(nextLine);
		// make sure nextLine starts with --
		if ((2 < strlen(nextLine))
				&&
			(nextLine[0] == '-')
			    &&
			(nextLine[1] == '-') ) {
			result = true;
			char *past = nextLine;
			while ((*past != '\0') && (!isspace(*past))) {
				++past;
			}
			if (*past != '\0') {
				*past = '\0';
				++past;
				past = puller.trim(past);
				past = strdup(past);
			}
			// prevent infinite loops caused by option files using
			//   optionFile option but allow a reasonable level of nesting
			if (0 == strcmp(nextLine, optionFile)) {
				++optFileDepth;
				if (optFileDepth < 10) {
					useOption(nextLine, past);
				}
				--optFileDepth;
			} else {
				useOption(nextLine, past);
			}
			nextLine = puller.getLine();
		} else {
			if (!result) {
				++lineCount;
				if (lineCount < 1000) {
					nextLine = puller.getLine();
				} else {
					nextLine = NULL;
				}
			}
		}
	}
	return result;
}




void OptionPass::getOptionsFromFile(const char *filePath) {
	LinePull puller(0x10000, filePath);
	char* nextLine = puller.getLine();
	while (nextLine != NULL) {
		nextLine = puller.trim(nextLine);
		// make sure nextLine starts with --
		if ((2 < strlen(nextLine))
				&&
			(nextLine[0] == '-')
			    &&
			(nextLine[1] == '-') ) {
			char *past = nextLine;
			while ((*past != '\0') && (!isspace(*past))) {
				++past;
			}
			if (*past != '\0') {
				*past = '\0';
				++past;
				past = puller.trim(past);
				past = strdup(past);
			}
			// prevent infinite loops caused by option files using
			//   optionFile option but allow a reasonable level of nesting
			if (0 == strcmp(nextLine, optionFile)) {
				++optFileDepth;
				if (optFileDepth < 10) {
					useOption(nextLine, past);
				}
				--optFileDepth;
			} else {
				useOption(nextLine, past);
			}
		}
		nextLine = puller.getLine();
	}
}

unsigned OptionPass::useOption(const char* optionType,
							 const char *optionValue) {
	// determine which option was specified and set appropriate values

	//  if optionType consumes optionValue than result will be incremented
	unsigned result = 0;
	if (strcmp(optionType, optionFile) == 0) {
		getOptionsFromFile(optionValue);
		++result;
	} else if (strcmp(optionType, pCutOpt) == 0) {
		gotOrthoSetting = true;
		eText = optionValue;
		maxLogE = getLog10(eText);
		++result;
	} else if (strcmp(optionType, piCutOpt) == 0) {
		gotOrthoSetting = true;
		identityCut = atof(optionValue);
		++result;
	} else if (strcmp(optionType, matchCutOpt) == 0) {
		gotOrthoSetting = true;
		matchCut = atof(optionValue);
		matchCut *= 0.01;
		++result;
	} else if (strcmp(optionType, maxWeightOpt) == 0) {
		gotOrthoSetting = true;
		maxWeight = atof(optionValue);
		++result;
	} else if (strcmp(optionType, linkOpt) == 0) {
		// -m 8 file was specified
		linkFile = optionValue;
		//  override any earlier bpo specification
		bpoInput = false;
		++result;
	} else if (strcmp(optionType, ggOpt) == 0) {
		ggFile = optionValue;
		++result;
	} else if (strcmp(optionType, workingOpt) == 0) {
		workingDir = optionValue;
		++result;
	} else if (strcmp(optionType, projectOpt) == 0) {
		projectName = optionValue;
		++result;
	} else if (strcmp(optionType, inflationOpt) == 0) {
		inflation = optionValue;
		++result;
	} else if (strcmp(optionType, cpuOpt) == 0) {
		cpuCount = optionValue;
		++result;
	} else if (strcmp(optionType, bOpt) == 0) {
		bBlast = optionValue;
		++result;
	} else if (strcmp(optionType, eOpt) == 0) {
		eBlast = optionValue;
		++result;
	} else if (strcmp(optionType, vOpt) == 0) {
		vBlast = optionValue;
		++result;
	} else if (strcmp(optionType, formatDbPath) == 0) {
		formatdb = optionValue;
		++result;
	} else if (strcmp(optionType, blastAllPath) == 0) {
		blastall = optionValue;
		++result;
	} else if (strcmp(optionType, mclPath) == 0) {
		mcl = optionValue;
		++result;
	} else if (strcmp(optionType, resultOpt) == 0) {
		setResultPath = true;
		resultFile = optionValue;
		++result;
	} else if (strcmp(optionType, setAll) == 0) {
		if (mixedCount == 0) {
			// first entry requires room for a single pointer
			mixedFastas = (const char **)(malloc(sizeof(char *)));
		} else {
			//  expand to create enough room to store next entry
			mixedFastas =
				(const char **)
					(realloc(mixedFastas, (mixedCount + 1) * sizeof(char *)));
		}
		mixedFastas[mixedCount] = optionValue;
		++mixedCount;
		++result;
	} else if (strcmp(optionType, addGenome) == 0) {
		if (singleCount == 0) {
			genomeFastas = (const char **) (malloc(sizeof(char *)));
		} else {
			genomeFastas =
				(const char **)
				   (realloc(genomeFastas, (singleCount + 1) * sizeof(char *)));
		}
		genomeFastas[singleCount] = optionValue;
		++singleCount;
		++result;
	} else if (strcmp(optionType, queryIndex) == 0) {
		indexSettings[QUERY_AT] = atoi(optionValue);
		++result;
	} else if (strcmp(optionType, subjectIndex) == 0) {
		indexSettings[SUBJECT_AT] = atoi(optionValue);
		++result;
	} else if (strcmp(optionType, eValueIndex) == 0) {
		indexSettings[E_VALUE_AT] = atoi(optionValue);
		++result;
	} else if (strcmp(optionType, identityIndex) == 0) {
		indexSettings[IDENTITY_AT] = atoi(optionValue);
		++result;
	} else if (strcmp(optionType, areaIndex) == 0) {
		indexSettings[LENGTH_AT] = atoi(optionValue);
		++result;
	} else if (strcmp(optionType, queryStartIndex) == 0) {
		indexSettings[QUERY_START] = atoi(optionValue);
		++result;
	} else if (strcmp(optionType, queryEndIndex) == 0) {
		indexSettings[QUERY_END] = atoi(optionValue);
		++result;
	} else if (strcmp(optionType, queryLengthIndex) == 0) {
		indexSettings[QUERY_LENGTH] = atoi(optionValue);
		++result;
	} else if (strcmp(optionType, subjectStartIndex) == 0) {
		indexSettings[SUBJECT_START] = atoi(optionValue);
		++result;
	} else if (strcmp(optionType, subjectEndIndex) == 0) {
		indexSettings[SUBJECT_END] = atoi(optionValue);
		++result;
	} else if (strcmp(optionType, subjectLengthIndex) == 0) {
		indexSettings[SUBJECT_LENGTH] = atoi(optionValue);
		++result;
	} else if (strcmp(optionType, mapColumnIndex) == 0) {
		indexSettings[MAP_COLUMN] = atoi(optionValue);
		++result;
	} else if (strcmp(optionType, tabSplit) == 0) {
		splitter = '\t';
	} else if (strcmp(optionType, bpoOpt) == 0) {
		bpoInput = true;
		// set defaults for bpo files
		splitter = ';';
		indexSettings[QUERY_AT] = 1;
		indexSettings[QUERY_LENGTH] = 2;
		indexSettings[SUBJECT_AT] = 3;
		indexSettings[SUBJECT_LENGTH] = 4;
	    indexSettings[E_VALUE_AT] = 5;
		indexSettings[IDENTITY_AT] = 6;
		indexSettings[MAP_COLUMN] = 7;
		linkFile = optionValue;
		++result;
	} else if (0 == strcmp(optionType, splitOpt)) {
		if (strlen(optionValue) == 1) {
			splitter = *optionValue;
		}
		++result;
	} else if (strcmp(optionType, orthoQuirkOpt) == 0) {
		gotOrthoSetting = true;
		useQuirks = true;
	} else if (strcmp(optionType, preBlastOpt) == 0) {
		blastPrep = true;
	} else if (strcmp(optionType, legacyOpt) == 0) {
		useLegacy = true;
	} else {
		cerr << "Invalid option " << optionType << endl;
	}
	return result;
}

char *OptionPass::getProjectFile(const char *suffix) {
	char *result = (char *) (malloc(projectFileFront + 1 + strlen(suffix)));
	memcpy(result, projectFileStart, projectFileFront);
	strcpy(result + projectFileFront, suffix);
	return result;
}

char *OptionPass::getMatrixPath() {
	if (matrixPath == NULL) {
		matrixPath = getProjectFile(".mtx");
	}
	return matrixPath;
}

char *OptionPass::getMclOutputPath() {
	return mclOutputPath;
}

GgParser * OptionPass::getGgParser() {
	return ggParser;
}

// used to begin process of saving options used
FILE* OptionPass::getOptionFile() {
	char *optFile = getProjectFile(".opt");
	FILE* optSave = fopen(optFile, "w");
	fprintf(optSave, "%s %s\n", workingOpt, workingDir);
	fprintf(optSave, "%s %s\n", projectOpt, projectName);
	return optSave;
}

bool checkFileDir(const char *filePath) {
	// remove file name from path
	int pathLen = strlen(filePath);
	// check to see if path is too short
	bool result = (1 < pathLen);
	if (result) {
		const char *endLook = filePath + pathLen - 1;
		// drop path separator at end of line
		if ((*endLook == '\\') || (*endLook == '/')) {
			--endLook;
			if ((*endLook == '\\') || (*endLook == '/')) {
				endLook = filePath;
			}
		}
		while ((filePath < endLook)
				  &&
			   ((*endLook != '\\') && (*endLook != '/'))) {
			--endLook;
		}
		result = (filePath < endLook);
		if (result) {
			pathLen = endLook - filePath;
			char *tmpPath = (char *)(malloc(pathLen + 1));
			memcpy(tmpPath, filePath, pathLen);
			tmpPath[pathLen] = '\0';
			struct stat statGet;
			if (stat(tmpPath, &statGet) == 0) {
				// tmpPath is at least a path
				DIR *dirHandle = opendir(tmpPath);
				if (dirHandle == NULL) {
					std::cerr << tmpPath << " is not a directory" << std::endl;
					result = false;
				} else {
					//  release handle on check directory
					closedir(dirHandle);
				}
			} else {
				result = false;
			}
			free(tmpPath);
		}
	}
	return result;
}

char *checkDir(const char *check) {
	char *result = NULL;
	if (check == NULL) {
		std::cerr << "No working directory setting" << std::endl;
	} else {
		int pathLen = strlen(check);
		char endChar = check[pathLen - 1];
		// drop path separator at end of line
		if ((endChar == '\\') || (endChar == '/')) {
			--pathLen;
		}
		result = (char *)malloc(pathLen + 1);
		memcpy(result, check, pathLen);
		result[pathLen] = '\0';
		// determine if check describes an existing directory
		struct stat statGet;
		if (stat(result, &statGet) == 0) {
			// check is at least a path
			DIR *dirHandle = opendir(check);
			if (dirHandle == NULL) {
				std::cerr << result << " is not a directory" << std::endl;
				result = NULL;
			} else {
				//  release handle on check directory
				closedir(dirHandle);
			}
		} else if (MKDIR != 0) {
			// unable to create a directory based on non-existent check path
			std::cerr << "Unable to create " << check << std::endl;
			result = NULL;
		}
	}
	return result;
}

const char *checkExist(const char *toCheck) {
	if (toCheck != NULL) {
		struct STAT64 statGet;
		int fileStatus = STAT64(toCheck, &statGet);
		if (fileStatus != 0) {
			printf("%d ", fileStatus);
			std::cerr << toCheck << " does not exist" << std::endl;
			toCheck = NULL;
		}
	}
	return toCheck;
}

bool checkFastas(unsigned count, const char **list) {
	bool result = true;
	if (0 < count) {
		for (unsigned i = 0; i < count; i++) {
			const char * nullCheck = checkExist(list[i]);
			result &= (nullCheck != NULL);
		}
	}
	return result;
}

FILE* OptionPass::storeFastaLists() {
	FILE *optSave = NULL;
	if ((0 < singleCount) || (0 < mixedCount)) {
		optSave = getOptionFile();
		for (unsigned i = 0; i < singleCount; i++) {
			fprintf(optSave, "%s %s\n", addGenome, genomeFastas[i]);
		}
		for (unsigned i = 0; i < mixedCount; i++) {
			fprintf(optSave, "%s %s\n", setAll, mixedFastas[i]);
		}
		ggFile = getProjectFile(".glg");
		idxFile = getProjectFile(".idx");
		ggParser =
			new GgParser(singleCount, genomeFastas,
						mixedCount, mixedFastas,
						ggFile, idxFile, allFasta);
	}
	return optSave;
}


bool OptionPass::createGgFromFastas(bool result) {
	if (result) {
		FILE* optSave = storeFastaLists();
		if (optSave != NULL) {
			fprintf(optSave, "%s\n", preBlastOpt);
			fclose(optSave);

		} else {
			result = false;
			std::cerr << "No fasta file specifications." << std::endl;
		}
	}
	return result;
}

void OptionPass::saveLinkDetails(FILE* optSave) {
	if (bpoInput) {
		fprintf(optSave, "%s %s\n", bpoOpt, linkFile);
		fprintf(optSave, "%s %d\n", queryLengthIndex,
				indexSettings[QUERY_LENGTH]);
		fprintf(optSave, "%s %d\n", subjectLengthIndex,
				indexSettings[SUBJECT_LENGTH]);
		fprintf(optSave, "%s %d\n", mapColumnIndex,
				indexSettings[MAP_COLUMN]);
	} else {
		fprintf(optSave, "%s %s\n", linkOpt, linkFile);
		fprintf(optSave, "%s %d\n", queryStartIndex,
					indexSettings[QUERY_START]);
		fprintf(optSave, "%s %d\n", queryEndIndex,
				indexSettings[QUERY_END]);
		fprintf(optSave, "%s %d\n", subjectStartIndex,
				indexSettings[SUBJECT_START]);
		fprintf(optSave, "%s %d\n", subjectEndIndex,
				indexSettings[SUBJECT_END]);
		fprintf(optSave, "%s %d\n", areaIndex, indexSettings[LENGTH_AT]);
	}
	fprintf(optSave, "%s %d\n", queryIndex, indexSettings[QUERY_AT]);
	fprintf(optSave, "%s %d\n", subjectIndex, indexSettings[SUBJECT_AT]);
	fprintf(optSave, "%s %d\n", eValueIndex, indexSettings[E_VALUE_AT]);
	fprintf(optSave, "%s %d\n", identityIndex, indexSettings[IDENTITY_AT]);
	if (splitter == '\t') {
		fprintf(optSave, "%s\n", tabSplit);
	} else {
		fprintf(optSave, "%s %c\n", splitOpt, splitter);
	}
}

void OptionPass::saveSettings(FILE* optSave) {
	fprintf(optSave, "%s %s\n", pCutOpt, eText);
	fprintf(optSave, "%s %f\n", piCutOpt, identityCut);
	fprintf(optSave, "%s %f\n", matchCutOpt, matchCut);
	fprintf(optSave, "%s %f\n", maxWeightOpt, maxWeight);
	fprintf(optSave, "%s %s\n", mclPath, mcl);
	fprintf(optSave, "%s %s\n", inflationOpt, inflation);
	fprintf(optSave, "%s %s\n", resultOpt, resultFile);
	if (useQuirks) {
		fprintf(optSave, "%s\n", orthoQuirkOpt);
	}
	fclose(optSave);
}


bool OptionPass::startBlastWorkflow(bool result) {
	if (result) {
		FILE* optSave = storeFastaLists();
		if (optSave != NULL) {
			if (useLegacy) {
				fprintf(optSave, "%s\n", legacyOpt);
			}
			fprintf(optSave, "%s %s\n", formatDbPath, formatdb);
			fprintf(optSave, "%s %s\n", blastAllPath, formatdb);
			fprintf(optSave, "%s %s\n", cpuOpt, cpuCount);
			fprintf(optSave, "%s %s\n", bOpt, bBlast);
			fprintf(optSave, "%s %s\n", vOpt, vBlast);
			fprintf(optSave, "%s %s\n", eOpt, eBlast);
			saveSettings(optSave);
		} else {
			result = false;
			std::cerr << "No fasta file specifications." << std::endl;
		}
	}
	return result;
}


bool OptionPass::areOk() {
	allFasta = NULL;
	idxFile = NULL;
	// make sure all specified protein files exist
	bool result = checkFastas(singleCount, genomeFastas);
	result &= checkFastas(mixedCount, mixedFastas);

	// make sure a directory has been specified for temporary files
	workingDir = checkDir(workingDir);
	result &= (workingDir != NULL);
	if (projectName == NULL) {
		result = false;
		std::cerr << "Missing prefix for temporary files" << std::endl;
	} else if (result) {
		// create projectFileFront to support temporary file names
		//   and default names for unspecified files
		int workingFront = strlen(workingDir);
		projectFileFront = strlen(projectName);
		projectFileStart =
			(char *)malloc(1 + workingFront + projectFileFront);
		memcpy(projectFileStart, workingDir, workingFront);
		projectFileStart[workingFront] = PATH_BREAK;
		++workingFront;
		memcpy(projectFileStart + workingFront, projectName,
				projectFileFront);
		projectFileFront += workingFront;
	}
	if (blastPrep) {
		// only work is to create a combined fasta file for blast
		result = createGgFromFastas(result);
	} else if (result) {
		// have working directory and project name

		// make sure a file path is available for the final output
		if (resultFile == NULL) {
			resultFile = getProjectFile(".end");
		} else {
			result = checkFileDir(resultFile);
		}

		// make sure flow can use an mcl executable
		if (executer->runWhich(mcl) != 0) {
			result = false;
			std::cerr << "Could not find " << mcl << std::endl;
		}

		// check to see if work flow is required to run blast
		if (linkFile == NULL) {
			if (executer->runWhich(formatdb) != 0) {
				result = false;
				std::cerr << "Could not find " << formatdb << std::endl;
			}
			if (executer->runWhich(blastall) != 0) {
				result = false;
				std::cerr << "Could not find " << blastall << std::endl;
			}
			if (result) {
				allFasta = getProjectFile(".faa");
			}
			result &= startBlastWorkflow(result);
		} else {
			// make sure link file exists
			linkFile = checkExist(linkFile);
			result = (linkFile != NULL);
			if (result) {
				// at this point either a gg file or a list of
				//     protein files is required
				if (ggFile == NULL) {
					FILE *saveOpt = storeFastaLists();
					result = (saveOpt != NULL);
					if (result) {
						saveLinkDetails(saveOpt);
						saveSettings(saveOpt);
					}
				} else {
					ggFile = checkExist(ggFile);
					if (ggFile != NULL) {
						FILE* optSave = getOptionFile();
						fprintf(optSave, "%s %s\n", ggOpt, ggFile);
						saveLinkDetails(optSave);
						saveSettings(optSave);
						ggParser =
							new GgParser(ggFile, getProjectFile(".idx"));
					}
				}
			}
		}
	}
	return result;
}

void OptionPass::doQuirkSplit(LinkParser *parser, Links** savers) {
	if (useQuirks) {
		Ortho classifier(parser, savers);
	} else {
		Normal classifier(parser, savers);
	}
}

bool OptionPass::classifyLinks(Links **bins) {
	// no classification work if run only created combined file for blast
	bool result = !blastPrep;
	if (result) {
		// check to see if blast run is required
		if (linkFile == NULL) {
			linkFile = getProjectFile(".out");
			char *allFasta = getProjectFile(".faa");
			if (useLegacy) {
				result =
					executer->runLegacy(allFasta, formatdb, blastall, eBlast,
										linkFile, cpuCount, bBlast, vBlast);
			} else {
				result =
					executer->runBlast(allFasta, formatdb, blastall, eBlast,
									   linkFile, cpuCount, bBlast, vBlast);
			}
		}
		time_t startTime;
		time (&startTime);
		if (result) {
			//  split pre mcl work required by settings
			if (0.0 < matchCut) {
				if ((NULL == ggParser->geneLengths) && (!bpoInput)) {
					std::cerr << "no gene lengths for match test" << std::endl;
				} else {
					//  all three filters apply
					if (bpoInput) {
						// create object for reading link file
						BpoAllParser chopper(ggParser, linkFile, splitter,
											 maxWeight, indexSettings);
						// object for filterings lines
						PickAll linker(&chopper, ggParser,
										maxLogE, identityCut, matchCut);
						// apply algorithm to filtered lines
						doQuirkSplit(&linker, bins);
					} else {
						// create object for reading link file
						CharSplitParser chopper(splitter, linkFile, maxWeight,
												indexSettings);
						// object for filtering lines
						PickAll linker(&chopper, ggParser,
										maxLogE, identityCut, matchCut);
						// apply desired algorithm to filtered lines
						doQuirkSplit(&linker, bins);
					}
				}
			} else {
				// make sure unused columns do not increase
				//  the minimum required column count.
				indexSettings[MAP_COLUMN] = 0;
				indexSettings[QUERY_LENGTH] = 0;
				indexSettings[SUBJECT_LENGTH] = 0;
				indexSettings[QUERY_START] = 0;
				indexSettings[QUERY_END] = 0;
				indexSettings[SUBJECT_START] = 0;
				indexSettings[SUBJECT_END] = 0;
				if (0.0 < identityCut) {
					if (bpoInput) {
						// get link parser
						BpoAllParser chopper(ggParser, linkFile, splitter,
											 maxWeight, indexSettings);
						// get line filter
						PickPiE linker(&chopper, ggParser,
										maxLogE, identityCut);
						// apply desired algorithm
						doQuirkSplit(&linker, bins);
					} else {
						// get link parser
						CharSplitParser chopper(splitter, linkFile, maxWeight,
											    indexSettings);
						// get line filter
						PickPiE linker(&chopper, ggParser,
										maxLogE, identityCut);
						// apply desired algorithm
						doQuirkSplit(&linker, bins);
					}
				} else {
					// remove last unused column
					indexSettings[IDENTITY_AT] = 0;
					//  No special bpo logic is required since
					//    e-value is just a simple column  read
					CharSplitParser chopper(splitter, linkFile, maxWeight,
											 indexSettings);
					// get weakest line filter
					PickE linker(&chopper, ggParser, maxLogE);
					// apply desired algorithm
					doQuirkSplit(&linker, bins);
				}
			}
		}
		time_t endTime;
		time (&endTime);
		double duration = difftime (endTime, startTime);
		printf (" %.2lf to classify blast hits\n", duration );
	}
	return result;
}

bool OptionPass::runMcl() {
	mclOutputPath = getProjectFile(".ocl");
	return (0 == executer->runMcl(mcl, inflation, matrixPath, mclOutputPath));
}

OptionPass::~OptionPass() {
	delete (executer);
}
