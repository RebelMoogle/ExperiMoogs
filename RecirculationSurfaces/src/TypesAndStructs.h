#pragma once
#include <memory>
#include <vector>
#include <tuple>
#include <Eigen/Core>
#include <ofPolyline.h>
//#include <vclibs/include/math/VecN.hh>
//#include <vclibs/include/math/ode.hh>

typedef Eigen::Matrix<double, 5, 1> Vector5;
//typedef VC::math::ode::Solution<double, VC::math::VecN<double, 3>> Solution3D;

typedef std::tuple<ofPolyline, std::vector<double>> PathlineTimes;

struct mogSolution {

private:
    // only allow const access.
    // positions
    std::vector< Eigen::Vector3d > m_states;
    // times
    std::vector< double > m_times;

public:

    mogSolution(){}

    mogSolution(const std::vector<Eigen::Vector3d>& states, const std::vector<double>& times)
        : m_states(states), m_times(times)
    {}

    size_t size() const
    {
        assert(m_states.size() == m_times.size() && "Element count of both positions and times should be the same!");
        return m_states.size();
    }

    // const access so vectors don't get changed.
    const std::vector<Eigen::Vector3d>& y() const
    {
        return m_states; // .at does bounds checking
    }

    const std::vector<double>& t() const
    {
        return m_times; // .at does bounds checking
    }

    void append(const Eigen::Vector3d &state, double time)
    {
        // this doesn't seem exception safe...
        m_states.push_back(state);
        m_times.push_back(time);
    }
};

typedef mogSolution Solution3D;

struct mogObserver {

private:
    // only allow const access.
    // positions
    Solution3D& m_results;
    const unsigned int m_MaxSteps;
    unsigned int stepCount = 0;

public:
    mogObserver(Solution3D& results, const unsigned int MaxSteps = UINT_MAX)
        : m_results(results), m_MaxSteps(MaxSteps)
    {}

    void operator()(const Eigen::Vector3d &integratedPosition, double integratedTime)
    {
        m_results.append(integratedPosition, integratedTime);

        // check steps.
        stepCount += 1;
        if (stepCount > m_MaxSteps) {
            stepCount = 0;
            throw std::exception("Integration failed, maximum steps reached. \n"); //TODO: print time as well.
        }
    }
};

struct mogDistanceResult {
    Vector5 MinPosition;
    Eigen::Vector3d DistanceVector;
    std::unique_ptr<mogSolution> myIntegration = nullptr; // has to use move semantic. Can be nullptr if not used.

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

struct mogEigenVectorValues {
    Vector5 atPosition;
    std::vector<Vector5> eigenVectors;
    Vector5 eigenValues;
};