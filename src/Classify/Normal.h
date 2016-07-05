/*
 * normal.h
 *
 *  Created on: Jun 21, 2010
 *      Author: mscott
 */

#ifndef NORMAL_H_
#define NORMAL_H_

#include "../LinkParser/LinkParser.h"
#include "../List/Links.h"

class Normal {
private:
	Links** savers;
	unsigned query;
	unsigned queryTaxon;
	unsigned subject;
	float nextE;
	void addLine(unsigned type);
public:
	Normal(LinkParser *parser, Links** savers);
	virtual ~Normal();
};

#endif /* NORMAL_H_ */
