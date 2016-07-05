//============================================================================
// Name        : FastOrtho.cpp
// Author      :
// Version     :
// Copyright   : Your copyright notice
// Description : Hello World in C++, Ansi-style
//============================================================================

#include <iostream>
#include "FileIO/LinePull.h"
#include "OptionPass.h"
#include "OrthoTop.h"
#include <ctype.h>

int main(int argc, const char** argv) {
	OptionPass options(argc, argv);
	if (options.areOk()) {
		OrthoTop master(&options);
	} else {
		std::cerr << "Invalid options" << std::endl;
	}
	return 0;
}
