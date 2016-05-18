#pragma once
#include "ofMain.h"
#include <../ThirdParty/Eigen/Core>
#include <../ThirdParty/vclibs/include/math/VecN.hh>
#include <../ThirdParty/vclibs/include/math/ode.hh>

typedef Eigen::Matrix<double, 5, 1> Vector5;
typedef VC::math::ode::Solution<double, VC::math::VecN<double, 3>> Solution3D;

struct mogCloudEigenVals{
	ofMesh surfaceMesh = ofMesh();

	std::vector<Vector5> FirstZeroEigenVector;

	std::vector<Vector5> SecondZeroEigenVector;

	std::vector<Vector5> EigenValues;
};

struct mogEigenVectorValues {
    std::vector<Vector5> eigenVectors;
    Vector5 atPosition;
    Vector5 eigenValues;
};

struct mogDistanceResult {
    Vector5 MinPosition;
    Eigen::Vector3d DistanceVector;
    std::unique_ptr<Solution3D> myIntegration = nullptr; // has to use move semantic. Can be nullptr if not used.

    mogDistanceResult()
    {};
    mogDistanceResult(Vector5 position, Eigen::Vector3d distanceVector) : MinPosition(position), DistanceVector(distanceVector)
    {};

    // compares distance. (works as a key, if you want to read it that way. ;) )
    bool operator<(const mogDistanceResult& other) const
    {
        return DistanceVector.norm() < other.DistanceVector.norm();
    }

};

struct mogSurfaceResult {
    Vector5 SurfacePosition;
    Eigen::Vector3d Normal;
    Eigen::Vector3d DistanceVector;
};