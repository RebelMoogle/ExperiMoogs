#pragma once
#include "ofMain.h"
#include <../ThirdParty/Eigen/Core>

typedef Eigen::Matrix<double, 5, 1> Vector5;

struct mogCloudEigenVals{
	ofMesh surfaceMesh = ofMesh();

	std::vector<Vector5> FirstZeroEigenVector;

	std::vector<Vector5> SecondZeroEigenVector;

	std::vector<Vector5> EigenValues;
};