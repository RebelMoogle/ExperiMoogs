#pragma once

#include "ofMain.h"

class ImportPly
{
public:
	ImportPly();
	~ImportPly();
	static ofMesh ImportPlyEigenCloud(std::string filePathName);

};

