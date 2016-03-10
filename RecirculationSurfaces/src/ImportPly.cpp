#include "ImportPly.h"

ImportPly::ImportPly()
{
}


ImportPly::~ImportPly()
{
}

// imports ascii csv of Eigen Cloud values
// should consist of position x y z, t, tau, eigen1 x, y, z, t, tau, eigen2 x, y, z, t, tau and eigenvalues 0-4 per line. 
// extra values will be ignored.
mogCloudEigenVals ImportPly::ImportEigenCloud(const std::string filePathName)
{
	ofBuffer csvFileBuffer;
	if(ofFile::doesFileExist(ofToDataPath(filePathName)))
	{
		csvFileBuffer = ofBufferFromFile(ofToDataPath(filePathName));
	}
	else
	{
		ofLog(ofLogLevel::OF_LOG_ERROR, "FILE DOES NOT EXIST, please try again.");
		throw(FileNotFound("File could not be found: " + ofToDataPath(filePathName)));
	}

	mogCloudEigenVals resultCloud;
	resultCloud.surfaceMesh.setMode(OF_PRIMITIVE_POINTS);

	std::string currentLine = csvFileBuffer.getFirstLine();

	// parse ply file. 
	// first line should read ply
	if(!ofIsStringInString(currentLine, "RebelEigenCSV"))
	{
		ofLogError() << "File does not contain RebelEigenCSV on line 1: " << currentLine;
		throw (InvalidHeader("File does not start with RebelEigenCSV on Line 1: " + currentLine));
	}
	else
	{
		ofLogNotice() << "reading eigen values csv file.";
	}

	unsigned lineNumber = 0;
	while(!csvFileBuffer.isLastLine())
	{
		currentLine = csvFileBuffer.getNextLine();
		++lineNumber;
		
		auto lineValues = ofSplitString(currentLine, ",", true, true);
		// check number of elements on each line for any missing or extra info.
		if(lineValues.size() < 20)
		{
			ofLogError() << "Too few values in line, cannot fill vertex, skipping line nr. " << lineNumber;
			continue;
		}
		
		// position vector
		resultCloud.surfaceMesh.addVertex(ofVec3f(ofToDouble(lineValues[0]), ofToDouble(lineValues[1]), ofToDouble(lineValues[2])));
		resultCloud.surfaceMesh.addTexCoord(ofVec2f(ofToDouble(lineValues[3]), ofToDouble(lineValues[4]))); // change to texcoords if necessary. 

		// eigenVector 1
		resultCloud.FirstZeroEigenVector.push_back((Vector5() << ofToDouble(lineValues[5]), ofToDouble(lineValues[6]), ofToDouble(lineValues[7]), ofToDouble(lineValues[8]), ofToDouble(lineValues[9])).finished());

		// eigenVector 2
		resultCloud.SecondZeroEigenVector.push_back((Vector5() << ofToDouble(lineValues[10]), ofToDouble(lineValues[11]), ofToDouble(lineValues[12]), ofToDouble(lineValues[13]), ofToDouble(lineValues[14])).finished());

		// eigenValues
		resultCloud.EigenValues.push_back((Vector5() << ofToDouble(lineValues[15]), ofToDouble(lineValues[16]), ofToDouble(lineValues[17]), ofToDouble(lineValues[18]), ofToDouble(lineValues[19])).finished());

	}
	
	return resultCloud;
}
