#include "ImportPly.h"

ImportPly::ImportPly()
{
}


ImportPly::~ImportPly()
{
}

ofMesh ImportPly::ImportPlyEigenCloud(const std::string filePathName)
{
	ofFile plyFile;
	if(plyFile.doesFileExist(ofToDataPath(filePathName)))
	{
		plyFile.open(ofToDataPath(filePathName));
	}
	else
	{
		ofLog(ofLogLevel::OF_LOG_ERROR, "FILE DOES NOT EXIST, please try again.");
		throw("NotFound");
	}
	
	auto plyBuffer = plyFile.readToBuffer();
	plyFile.close(); // fine to close here, right?

	std::string currentLine = plyBuffer.getFirstLine();

	// parse ply file. 
	// first line should read ply
	if(currentLine.compare("ply") != 0)
	{
		ofLogError() << "Invalid header on line 1: " << currentLine;
		throw ("InvalidHeader");
	}
	else
	{
		ofLogNotice() << "reading ply file.";
	}

	bool isHeader = true;
	while(!plyBuffer.isLastLine())
	{
		currentLine = plyBuffer.getNextLine();
		if(!ofIsStringInString(currentLine, "format") && ofIsStringInString(currentLine, "ascii"))
		{
			// ascii, all god
			ofLogNotice() << currentLine;
			continue;
		}
		else
		{
			ofLogError() << "Importer currently only supports ASCII format, invalid Format: " << currentLine;
			throw("FormatNotSupported");
		}

		if(currentLine.compare())
		// next line should read file type. 
		// there might be arbitrary comments
		// then properties, 
		// header ends, data begins. 
		// check number of elements on each line for any missing or extra info.
	}




	return ofMesh();
}
