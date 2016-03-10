#pragma once

#include "ofMain.h"
#include "TypesAndStructs.h"

class FileNotFound : public std::runtime_error
{
public:
	FileNotFound(std::string message) : std::runtime_error(message) {}
	FileNotFound() : std::runtime_error("File not found") {}
};

class InvalidHeader : public std::runtime_error
{
public:
	InvalidHeader(std::string message) : std::runtime_error(message) {}
	InvalidHeader() : std::runtime_error("Invalid Header.") {}
};

class ImportPly
{
public:
	ImportPly();
	~ImportPly();
	static mogCloudEigenVals ImportEigenCloud(std::string filePathName);

};

