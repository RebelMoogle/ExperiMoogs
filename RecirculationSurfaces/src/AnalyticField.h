#pragma once
#include <functional>
#include <vector>
#include "TypesAndStructs.h"
//#include <vclibs/include/math/rk43.hh>
//#include <vclibs/include/math/VecN.hh>

class AnalyticField
{
public:
	AnalyticField(std::string fieldName, std::function<Eigen::Vector3d(const Eigen::Vector3d& x, Eigen::Vector3d& dxdt, const double t)> analyticFormula);
	~AnalyticField();

    inline void operator()(const Eigen::Vector3d& position, Eigen::Vector3d& resultPosition, const double t) const
    {
        AnalyticFormula(position, resultPosition, t);
    }

    inline Eigen::Vector3d GetDataAt(const Eigen::Vector3d& position, double t) const
    {
        Eigen::Vector3d result;
        (*this)(position, result, t);
        return result;
    }

	// Calculates Pathline starting from WorldPosition (includes time)
	// Can try to compute outside of bounding box if IgnoreBounds is true.
	//TODO make FRunnable to perform in threads
	// FORWARD ONLY
	void ComputePathLineAt(const Eigen::Vector4d& WorldPosition, std::vector<Eigen::Vector3d>& PathLine, const bool IgnoreBounds = false, const double StepSize = 0.1, const float IntegrationTime = FLT_MAX, const unsigned int maxSteps = 100000);

    const Eigen::Matrix<double, 3, 4> ComputeFlowGradientTime(const Eigen::Vector4d& position, const double IntegrationTime, const double StepSize = 0.001, const double CellSize = 0.0001) const;

    // adaptive RK43 integrator using boost odeint
    Eigen::Vector3d IntegrateOverTime(const Eigen::Vector4d& WorldPosition, const double IntegrationTime, const double StepSize = 0.01, Solution3D* Solution = nullptr, double absoluteError = 1.0e-10, double relativeError = 1.0e-6) const;
    
    inline Eigen::Vector3d IntegrateOverTime(const Vector5& WorldPosition, const double StepSize = 0.01, Solution3D* Solution = nullptr, double absoluteError = 1.0e-10, double relativeError = 1.0e-6) const
    {
        return IntegrateOverTime((Eigen::Vector4d() << WorldPosition(0), WorldPosition(1), WorldPosition(2), WorldPosition(3)).finished(), WorldPosition(4), StepSize, Solution, absoluteError, relativeError);
    }

	virtual inline const std::string GetName() const
	{
		return FieldName;
	}

private:
	const std::function<Eigen::Vector3d(const Eigen::Vector3d& x, Eigen::Vector3d& dxdt, const double t)> AnalyticFormula;

	const std::string FieldName;

};

