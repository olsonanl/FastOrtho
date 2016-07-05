/*
 * Executer.h
 *
 *  Created on: Jun 4, 2010
 *      Author: mscott
 */

#ifndef EXECUTER_H_
#define EXECUTER_H_

class Executer {
	char *commandBuffer;
	int commandLen;
public:
	Executer();
	void startCommand(const char *intro);
	void appendToCommand(const char *nextPart);
	int execute();
	int addLast(const char *final);
	bool runBlast(const char *allFasta, const char *formatdb,
				 const char *blastall, const char *eText,
				 const char *linkFile, const char *cpuCount,
				 const char *bBlast, const char *vBlast);
	bool runLegacy(const char *allFasta, const char *formatdb,
				    const char *blastall, const char *eText,
				    const char *linkFile, const char *cpuCount,
				    const char *bBlast, const char *vBlast);
	bool runMcl(const char *mcl, const char *inflation, const char *matrixFile,
				const char *clusterFile);

	int runWhich(const char *toCheck);
	virtual ~Executer();
};

#endif /* EXECUTER_H_ */
