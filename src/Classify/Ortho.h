/*
 * Ortho.h
 *
 *  Created on: Jun 21, 2010
 *      Author: mscott
 */

#ifndef ORTHO_H_
#define ORTHO_H_
#include "../LinkParser/LinkParser.h"
#include "../List/Links.h"
#include "../List/Blocks.h"


class Ortho {
private:
	Links** savers;

	float nextE;
	unsigned taxonCount;

	unsigned intraCount;
	Blocks *intraSave;

	unsigned query;
	unsigned queryTaxon;
	unsigned subject;
	unsigned subjectTaxon;

	unsigned bitBytes;
	unsigned char *saved;

	float *minInters;
	float minInter;

	void addLine(unsigned type, float eValue);
	void setNewQuery(LinkParser *parser);
	void storeIntra();

public:
	Ortho(LinkParser *parser, Links** savers);
	virtual ~Ortho();
};

#endif /* ORTHO_H_ */
