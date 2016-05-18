#pragma once
#include <functional>
#include <vector>
#include "TypesAndStructs.h"
#include <../ThirdParty/vclibs/include/math/rk43.hh>
#include <../ThirdParty/vclibs/include/math/VecN.hh>

class AnalyticField
{
public:
	AnalyticField(std::string fieldName, std::function<Eigen::Vector3d(const Eigen::Vector4d position)> analyticFormula);
	~AnalyticField();

	Eigen::Vector3d GetDataAt(Eigen::Vector4d position) const;

	// Calculates Pathline starting from WorldPosition (includes time)
	// Can try to compute outside of bounding box if IgnoreBounds is true.
	//TODO make FRunnable to perform in threads
	// FORWARD ONLY
	void ComputePathLineAt(const Eigen::Vector4d& WorldPosition, std::vector<Eigen::Vector3d>& PathLine, const bool IgnoreBounds = false, const double StepSize = 0.1, const float IntegrationTime = FLT_MAX, const unsigned int maxSteps = 100000);

    const Eigen::Matrix<double, 3, 4> ComputeFlowGradientTime(const Eigen::Vector4d position, const float IntegrationTime, const double StepSize = 0.001, const double CellSize = 0.0001) const;

	// adaptive RK43 integrator from vclibs
	Eigen::Vector3d IntegrateOverTimeVCLibsRK43(const Eigen::Vector4d& WorldPosition, const double IntegrationTime, const double StepSize = 0.01, VC::math::ode::Solution<double, VC::math::VecN<double, 3>>* Solution = nullptr, int MaximumSteps = 0) const;
	inline Eigen::Vector3d IntegrateOverTimeVCLibsRK43(const Vector5& WorldPosition, const double StepSize = 0.01, VC::math::ode::Solution<double, VC::math::VecN<double, 3>>* Solution = nullptr, int MaximumSteps = 0) const
	{
		return IntegrateOverTimeVCLibsRK43((Eigen::Vector4d() << WorldPosition(0), WorldPosition(1), WorldPosition(2), WorldPosition(3)).finished(), WorldPosition(4), StepSize, Solution, MaximumSteps);
	}

	virtual inline const std::string GetName() const
	{
		return FieldName;
	}

private:
	const std::function<Eigen::Vector3d(const Eigen::Vector4d position)> AnalyticFormula;

	const std::string FieldName;

};

