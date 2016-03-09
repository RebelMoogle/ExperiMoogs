#pragma once

#include "ofMain.h"
#include "TypesAndStructs.h"

class ImportPly
{
public:
	ImportPly();
	~ImportPly();
	static mogCloudEigenVals ImportEigenCloud(std::string filePathName);

};

