/*
 * Executer.cpp
 *
 *  Created on: Jun 4, 2010
 *      Author: mscott
 */

#include "Executer.h"
#include <stdlib.h>
#include <string.h>
#include <iostream>


Executer::Executer() {
	commandBuffer = (char *)(malloc(0x1000 * sizeof(char)));
	commandLen = 0;
}

void Executer::startCommand(const char *intro) {
	commandLen = strlen(intro);
	memcpy(commandBuffer, intro, commandLen);
}
void Executer::appendToCommand(const char *append) {
	int partSize = strlen(append);
	memcpy(commandBuffer + commandLen, append, partSize);
	commandLen += partSize;
}

int Executer::execute() {
	commandBuffer[commandLen] = '\0';
	std::cout << commandBuffer << std::endl;
	return (system(commandBuffer));
}

int Executer::addLast(const char *tail) {
	strcpy(commandBuffer + commandLen, tail);
	std::cout << commandBuffer << std::endl;
	return (system(commandBuffer));
}

int Executer::runWhich(const char *toFind) {
	int result;
	startCommand("which ");
	if (memcmp(toFind, "C:", 2) != 0) {
		result = addLast(toFind);
	} else {
		appendToCommand("/cygdrive/c");
		int len = strlen(toFind);
		for (int i = 2; i < len; i++) {
			char addChar = toFind[i];
			if (addChar == '\\') {
				addChar = '/';
			}
			commandBuffer[commandLen] = addChar;
			++commandLen;
		}
		result = execute();
	}
	return result;
}

bool Executer::runBlast(const char *allFasta, const char *formatdb,
						const char *blastall, const char *eText,
						const char *linkFile, const char *cpuCount,
						const char *bBlast, const char *vBlast) {
	startCommand(formatdb);
	appendToCommand(" -dbtype prot -in ");
	bool result = (addLast(allFasta) == 0);
	if (result) {
		startCommand(blastall);
		appendToCommand(" -db ");
		appendToCommand(allFasta);
		appendToCommand(" -query ");
		appendToCommand(allFasta);
		appendToCommand(" -evalue ");
		appendToCommand(eText);
		appendToCommand(" -num_threads ");
		appendToCommand(cpuCount);
		appendToCommand(" -num_descriptions ");
		appendToCommand(vBlast);
		appendToCommand(" -num_alignments ");
		appendToCommand(bBlast);
		appendToCommand(" -out ");
		appendToCommand(linkFile);
		result = (addLast(" -outfmt 6") == 0);
	}
	return result;
}

bool Executer::runLegacy(const char *allFasta, const char *formatdb,
						 const char *blastall, const char *eText,
						 const char *linkFile, const char *cpuCount,
						 const char *bBlast, const char *vBlast) {
	startCommand(formatdb);
	appendToCommand(" -i ");
	appendToCommand(allFasta);
	bool result = (addLast(" -p t") == 0);
	if (result) {
		startCommand(blastall);
		appendToCommand(" -p blastp -i ");
		appendToCommand(allFasta);
		appendToCommand(" -d ");
		appendToCommand(allFasta);
		appendToCommand(" -e ");
		appendToCommand(eText);
		appendToCommand(" -o ");
		appendToCommand(linkFile);
		appendToCommand(" -m 8 -a ");
		appendToCommand(cpuCount);
		appendToCommand(" -v ");
		appendToCommand(vBlast);
		appendToCommand(" -b ");
		result = (addLast(bBlast) == 0);
	}
	return result;
}


bool Executer::runMcl(const char *mcl, const char *inflation,
					  const char *matrixFile, const char *clusterFile) {
	startCommand(mcl);
	appendToCommand(" ");
	appendToCommand(matrixFile);
	appendToCommand(" -I ");
	appendToCommand(inflation);
	appendToCommand(" -o ");
	return (addLast(clusterFile));
}

Executer::~Executer() {
	free(commandBuffer);
}
