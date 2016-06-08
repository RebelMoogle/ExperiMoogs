#include "AnalyticField.h"
#include "ofPolyline.h"

#include <boost/numeric/odeint.hpp>
#include <boost/numeric/odeint/external/eigen/eigen.hpp>

AnalyticField::AnalyticField(std::string fieldName, std::function<Eigen::Vector3d(const Eigen::Vector3d& x, Eigen::Vector3d& dxdt, const double t)> analyticFormula):FieldName(fieldName), AnalyticFormula(analyticFormula)
{

}


AnalyticField::~AnalyticField()
{
}

void AnalyticField::ComputePathLineAt(const Eigen::Vector4d & WorldPosition, ofPolyline& PathLine, const bool IgnoreBounds, const double StepSize, const float IntegrationTime, const unsigned int maxSteps)
{
    std::clog << "Computing PathLine at " << WorldPosition << "\n";

    Eigen::Vector3d CurPos3;
    unsigned int CurrentStep = 0;

    // Add starting / and (hopefully) final position: 
    PathLine.addVertex(WorldPosition[0], WorldPosition[1], WorldPosition[2]);

    Solution3D Solution;
    IntegrateOverTime(WorldPosition, IntegrationTime, StepSize, &Solution);

    // add positions of solution ofpolyline (of type)
    for(auto position : Solution.y()) 
    {
        PathLine.addVertex(position(0), position(1), position(2));
    }

}

const Eigen::Matrix<double, 3, 4> AnalyticField::ComputeFlowGradientTime(const Eigen::Vector4d& position, const double IntegrationTime, const double StepSize, const double CellSize) const
{
    // TODO: integrate over timestep ( pass time for integration to function!)
    Eigen::Vector3d XHi = IntegrateOverTime(position + Eigen::Vector4d( CellSize, 0., 0., 0.), IntegrationTime, StepSize);
    Eigen::Vector3d XLo = IntegrateOverTime(position + Eigen::Vector4d(-CellSize, 0., 0., 0.), IntegrationTime, StepSize);
    Eigen::Vector3d YHi = IntegrateOverTime(position + Eigen::Vector4d(0.,  CellSize, 0., 0.), IntegrationTime, StepSize);
    Eigen::Vector3d YLo = IntegrateOverTime(position + Eigen::Vector4d(0., -CellSize, 0., 0.), IntegrationTime, StepSize);
    Eigen::Vector3d ZHi = IntegrateOverTime(position + Eigen::Vector4d(0., 0.,  CellSize, 0.), IntegrationTime, StepSize);
    Eigen::Vector3d ZLo = IntegrateOverTime(position + Eigen::Vector4d(0., 0., -CellSize, 0.), IntegrationTime, StepSize);
    Eigen::Vector3d THi = IntegrateOverTime(position + Eigen::Vector4d(0., 0., 0.,  CellSize), IntegrationTime, StepSize);
    Eigen::Vector3d TLo = IntegrateOverTime(position + Eigen::Vector4d(0., 0., 0., -CellSize), IntegrationTime, StepSize);


    Eigen::Vector3d dx = (XHi - XLo) / (2.0*CellSize);
    Eigen::Vector3d dy = (YHi - YLo) / (2.0*CellSize);
    Eigen::Vector3d dz = (ZHi - ZLo) / (2.0*CellSize);
    Eigen::Vector3d dt = (THi - TLo) / (2.0*CellSize);

    Eigen::Matrix<double, 3, 4> ResultGradient;
    ResultGradient << dx, dy, dz, dt;
    return ResultGradient;
}

Eigen::Vector3d AnalyticField::IntegrateOverTime(const Eigen::Vector4d& WorldPosition, const double IntegrationTime, const double StepSize, Solution3D* Solution, const int MaxSteps) const
{
    namespace boost_ode = boost::numeric::odeint;
    //integrate_adaptive only returns number of steps. needs observer to return something. 
    // if we are only interested in the last value, use custom observer that only saves the last value. 
     // default = 500



    if (Solution == nullptr)         
    {
        Eigen::Vector3d finalPosition;
        int stepCount = 0;
        // function call takes: Integrator, integration function, start position, start time, end time, step size, observer function (save result)
        //returns number of steps as size_t.
        try {
            boost_ode::integrate<double>(*this, Eigen::Vector3d(WorldPosition[0], WorldPosition[1], WorldPosition[2]), WorldPosition[3], WorldPosition[3] + IntegrationTime, StepSize,
                [&finalPosition, &stepCount, &MaxSteps](const Eigen::Vector3d &integratedPosition, double integratedTime)
            {
                finalPosition = integratedPosition;
                stepCount += 1;
                if (stepCount > MaxSteps) {
                    throw std::exception("Integration failed: Maximum step count reached. \n"); // TODO: print current time as well.
                }
            });
        }
        catch (const std::exception& e) {
            std::clog << e.what();
        }

        // uses Lambda function to save

        return finalPosition;
    }
    else         
    {
        try {
            //returns number of steps as size_t. (could use for debugging info)
            boost_ode::integrate<double>(*this, Eigen::Vector3d(WorldPosition[0], WorldPosition[1], WorldPosition[2]), WorldPosition[3], WorldPosition[3] + IntegrationTime, StepSize, mogObserver(*Solution, MaxSteps));
        }
        catch (const std::exception& e) {
           std::clog << e.what();
        }


        return Solution->y().back();
    }
    

    return Eigen::Vector3d();
}