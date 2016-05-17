#include "AnalyticField.h"


class Evaluator {
public:
	typedef double real_t;
	typedef VC::math::VecN<double, 3> vec_t;

	Evaluator() :FlowData(nullptr) {};
	Evaluator(const AnalyticField* flowData) : FlowData(flowData) {};

	bool refine(vec_t&, vec_t&) {
		return false; // we dont want / need to refine this vector field
	}

	void initialize_vector(vec_t& _x, const vec_t& _y) {}

	void initialize(real_t, const vec_t&,
		VC::math::ode::RK43<Evaluator>*) {
	}

	void dy(real_t _t, const vec_t& _y, vec_t& _dy)
	{
		assert(flowData != nullptr);
		Eigen::Vector3d Velocity = FlowData->GetDataAt((Eigen::Vector4d() << _y[0], _y[1], _y[2], _t).finished());
		_dy = vec_t(Velocity(0), Velocity(1), Velocity(2));
	}

	void output(real_t _t, const vec_t& _y, const vec_t& _dy)
	{
		// log each step.?
		// action to perform / output each step
	}

	bool is_inside(real_t, const vec_t&) {
		return true;
	}

private:
	const AnalyticField* FlowData;

};

AnalyticField::AnalyticField(std::string fieldName, std::function<Eigen::Vector3d(const Eigen::Vector4d position)> analyticFormula):FieldName(fieldName), AnalyticFormula(analyticFormula)
{

}


AnalyticField::~AnalyticField()
{
}

Eigen::Vector3d AnalyticField::GetDataAt(Eigen::Vector4d position) const
{
    return AnalyticFormula(position);
}

void AnalyticField::ComputePathLineAt(const Eigen::Vector4d & WorldPosition, std::vector<Eigen::Vector3d>& PathLine, const bool IgnoreBounds, const double StepSize, const float IntegrationTime, const unsigned int maxSteps)
{
    float EndTime = WorldPosition.w() + IntegrationTime;

    std::clog << "Computing PathLine at " << WorldPosition << "\n";

    Eigen::Vector4d CurrentPosition = WorldPosition;
    Eigen::Vector3d FinalDestination = WorldPosition;
    Eigen::Vector3d CurPos3;
    unsigned int CurrentStep = 0;

    // Add starting / and (hopefully) final position: 
    PathLine.push_back(FinalDestination);

    VC::math::ode::Solution<double, VC::math::VecN<double, 3>> Solution;
    IntegrateOverTimeVCLibsRK43(CurrentPosition, IntegrationTime, StepSize, &Solution);

    // add positions of solution to pathline array
    for_each(Solution.y.begin(), Solution.y.end(), [&](VC::math::VecN<double, 3> currentVector)
    {
        PathLine.push_back(Eigen::Vector3d(currentVector.data())); // create Eigen::Vector3d from vclibs vector
    });
}

Eigen::Vector3d AnalyticField::IntegrateOverTimeVCLibsRK43(const Eigen::Vector4d & WorldPosition, const double IntegrationTime, const double StepSize, VC::math::ode::Solution<double, VC::math::VecN<double, 3>>* Solution, int MaximumSteps) const
{
    VC::math::ode::RK43<Evaluator> rk43;
    VC::math::ode::Solution<double, VC::math::VecN<double, 3>>* Sol;
    VC::math::ode::EvalState state;
    Evaluator eval = Evaluator(this);
    rk43.options.hmax = 0.1;
    rk43.options.rsmin = 1e-24;

    if (Solution) {

        Sol = Solution;
    }
    else {
        Sol = new VC::math::ode::Solution<double, VC::math::VecN<double, 3>>();
    }


    if (MaximumSteps) {
        state = rk43.integrate(&eval, VC::math::VecN<double, 3>(WorldPosition(0), WorldPosition(1), WorldPosition(2)), WorldPosition(3), WorldPosition(3) + IntegrationTime, Sol, false, MaximumSteps);
    }
    else {
        state = rk43.integrate(&eval, VC::math::VecN<double, 3>(WorldPosition(0), WorldPosition(1), WorldPosition(2)), WorldPosition(3), WorldPosition(3) + IntegrationTime, Sol);
    }

    if (state && state != VC::math::ode::Success) {
        //UE_LOG(GeneralFC, Error, TEXT("Integration failed: %s"), ANSI_TO_TCHAR(VC::math::ode::StateStr[state]));

        return Eigen::Vector3d(rk43.y().data());
    }

    /*// output solution
    if (argc==0 || (argc>1 && argv[1][0]=='s'))
    cout << sol << endl;
    */
    if (!Solution) {
        delete Sol;
    }



    return (Eigen::Vector3d(rk43.y().data()));
}
