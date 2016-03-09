#pragma once

#include <Eigen/Core>

typedef Eigen::Matrix<double, 5, 1> Vector5;

struct mogCloudEigenVals
{
	ofMesh surfaceMesh;

	std::vector<Vector5> FirstZeroEigenVector;

	std::vector<Vector5> SecondZeroEigenVector;

	std::vector<Vector5> EigenValues;
};
