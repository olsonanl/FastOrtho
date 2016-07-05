/*
 * OrthoTop.h
 *
 *  Created on: Apr 8, 2010
 *      Author: mscott
 */

#ifndef ORTHOTOP_H_
#define ORTHOTOP_H_
#include "OptionPass.h"
#include "LinkParser/LinkParser.h"
#include "List/LinkFinder.h"
#include <vector>



class OrthoTop {
private:
	GgParser *ggParser;

	LinkFinder minCheck;
	LinkFinder maxCheck;
	LinkFinder toFind;

	double **interMeans;
	unsigned taxonCount;

	Links* intraTops;
	Links* bestInters;
	Links* otherInters;

	void setInterMeans(unsigned taxonCount, unsigned *getParent);
	void expandInters();
	void normalizeInters(unsigned *gParent);
	bool doHubCheck(unsigned low, std::vector<unsigned> highClones);
	bool doHubCheck(std::vector<unsigned> lowClones, unsigned high);
	bool viewHubCheck(std::vector<unsigned> lowClones, unsigned high);
	bool doHubCheck(std::vector<unsigned> lowClones,
			std::vector<unsigned> highClones);
	bool viewHubCheck(std::vector<unsigned> lowClones,
			std::vector<unsigned> highClones);
	void adjustInters(unsigned *getParent);
	void writeMclInput(char *matrixPath);


public:
	void doOrthoConvert(const char *mclOutput, const char *idxPath,
						const char *readable);
	OrthoTop(OptionPass *options);
	virtual ~OrthoTop();
};

#endif /* ORTHOTOP_H_ */
