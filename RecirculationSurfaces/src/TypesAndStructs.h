#pragma once
#include <memory>
#include <vector>
#include <Eigen/Core>
//#include <vclibs/include/math/VecN.hh>
//#include <vclibs/include/math/ode.hh>

typedef Eigen::Matrix<double, 5, 1> Vector5;
//typedef VC::math::ode::Solution<double, VC::math::VecN<double, 3>> Solution3D;
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

    int size() const
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

public:
    mogObserver(Solution3D& results)
        : m_results(results)
    {}

    void operator()(const Eigen::Vector3d &integratedPosition, double integratedTime)
    {
        m_results.append(integratedPosition, integratedTime);
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