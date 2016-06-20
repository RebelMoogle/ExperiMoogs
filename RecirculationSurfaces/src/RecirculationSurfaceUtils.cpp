#include "RecirculationSurfaceUtils.h"
#include <memory>
#include <Eigen/EigenValues>
#include <iostream>

RecirculationSurfaceUtils::~RecirculationSurfaceUtils()
{}

const mogDistanceResult RecirculationSurfaceUtils::FindSurfacePositionFrom(Vector5 StartPosition, double StepSize, double MinimumStep, const int maximumIterations)
{

    if (MinimumStep == 0.0) {
        MinimumStep = StepSize * 1.0e-2;
    }
    int numIterations = 0;
    mogDistanceResult result;
    result.MinPosition = StartPosition;
    Eigen::Vector3d integratedPosition = FlowData.IntegrateOverTime(StartPosition, StepSize);
    result.DistanceVector = integratedPosition - Eigen::Vector3d(StartPosition[0], StartPosition[1], StartPosition[2]);
    double Distance = result.DistanceVector.norm();

    if (Distance < parameters.MaxError) {
        //already close enough.
        return result;
    }


    //UE_LOG(GeneralFC, Log, TEXT("Starting to find closest position on surface, with a starting StepSize / Error of %f and starting Distance of %f"), StepSize, OldDistance);
    do {
        StartPosition += ComputeGradientDescentToSurface(StartPosition, integratedPosition, StepSize) * StepSize;
        integratedPosition = FlowData.IntegrateOverTime(StartPosition, StepSize);
        Eigen::Vector3d NewDistanceVector = integratedPosition - Eigen::Vector3d(StartPosition[0], StartPosition[1], StartPosition[2]);
        double NewDistance = NewDistanceVector.norm();

        if (NewDistance < parameters.MaxError) {
            // distance is good enough. (terminate early)
            result.MinPosition = StartPosition;
            result.DistanceVector = NewDistanceVector;
            Distance = NewDistance;
            return result;
        }
        else if ((NewDistance / StartPosition[4]) < (Distance / result.MinPosition[4])) {
            // got smaller, continue until it doesn't get smaller.
            result.MinPosition = StartPosition;
            result.DistanceVector = NewDistanceVector;
            Distance = NewDistance;
        }
        else if (NewDistance > parameters.MaxError) {
            //Error too large, decrease stepsize! //maybe, if we're lucky, it will help. Not likely though.
            StepSize = StepSize*0.1;
            StartPosition = result.MinPosition;
        }

    } while (StepSize > MinimumStep && maximumIterations >= numIterations++);
    

    return result; // OldPosition will have the best one we found so far.
}

const std::vector<mogDistanceResult>& RecirculationSurfaceUtils::GetDistances()
{
    if (MinimumFutureDistances.size() > 0) {
        MinimumDistances.clear();
        for (unsigned int i = 0; i < MinimumFutureDistances.size(); ++i) {
            MinimumDistances.push_back(MinimumFutureDistances[i].get());
        }
        MinimumFutureDistances.clear();
    }
    else if (MinimumDistances.size() == 0) {
        // no calculation started
        throw std::exception("Please start Distance Calculation first!");
    }

    return MinimumDistances;
}

const Vector5 RecirculationSurfaceUtils::ComputeGradientDescentToSurface(const Vector5 StartVector, const Eigen::Vector3d IntegratedPosition, const double StepSize) const
{
    //Flowmap at position - position
    const Eigen::Vector3d Distance = IntegratedPosition - Eigen::Vector3d(StartVector(0), StartVector(1), StartVector(2));
    const Eigen::Matrix<double, 3, 5> SurfaceMatrixM = CreateSurfaceMatrixM(StartVector, IntegratedPosition, StepSize);

    return (-SurfaceMatrixM.transpose() * Distance);
}

const Eigen::Matrix<double, 3, 5> RecirculationSurfaceUtils::CreateSurfaceMatrixM(const Vector5 SurfacePosition, const Eigen::Vector3d IntegratedPosition, double StepSize) const
{
    Eigen::Matrix<double, 3, 5>  resultMatrix;

    // 5th dimension is integration time at surface position.
    if (SurfacePosition[4] == 0.0f) {
        auto JacobianT = FlowData.ComputeFlowGradientTime(Eigen::Vector4d(SurfacePosition[0], SurfacePosition[1], SurfacePosition[2], SurfacePosition[3]), SurfacePosition[4], StepSize);
        Eigen::Matrix<double, 3, 3> Jacobian;
        Eigen::Vector3d curFlow = JacobianT.col(3);
        Jacobian << JacobianT.col(0), JacobianT.col(1), JacobianT.col(2);


        Eigen::Vector3d curVD2 = Jacobian * curFlow;

        resultMatrix << Jacobian(0, 0), Jacobian(0, 1), Jacobian(0, 2), curFlow(0, 3), curVD2(0),
            Jacobian(1, 0), Jacobian(1, 1), Jacobian(1, 2), curFlow(1, 3), curVD2(1),
            Jacobian(2, 0), Jacobian(2, 1), Jacobian(2, 2), curFlow(2, 3), curVD2(2);

        return resultMatrix;
    }


    Eigen::Matrix<double, 3, 3> IdentityMatrix;
    IdentityMatrix << 1, 0, 0,
        0, 1, 0,
        0, 0, 1;

    const Eigen::Matrix<double, 3, 4> GradientGT = FlowData.ComputeFlowGradientTime(Eigen::Vector4d(SurfacePosition[0], SurfacePosition[1], SurfacePosition[2], SurfacePosition[3]), SurfacePosition[4], StepSize);

    const Eigen::Matrix<double, 3, 3> GradientG = (Eigen::Matrix<double, 3, 3>() << GradientGT.col(0), GradientGT.col(1), GradientGT.col(2)).finished() - IdentityMatrix;

    //const Eigen::Matrix<double, 3, 3> GradientG = FlowData->ComputeFlowGradient(FVector4(SurfacePosition[0], SurfacePosition[1], SurfacePosition[2], SurfacePosition[3]), SurfacePosition[4], StepSize) - IdentityMatrix;

    const Eigen::Vector3d FirstPosEigen = (Eigen::Vector3d() << SurfacePosition[0], SurfacePosition[1], SurfacePosition[2]).finished();
    const Eigen::Vector3d VelocityDAtEnd = FlowData.GetDataAt(IntegratedPosition, SurfacePosition[3] + SurfacePosition[4]) - (IntegratedPosition - FirstPosEigen) / SurfacePosition[4];

    //[row][column] - for both Eigen and Unreal
    resultMatrix << GradientG(0, 0), GradientG(0, 1), GradientG(0, 2), GradientGT(0, 3), VelocityDAtEnd(0),
        GradientG(1, 0), GradientG(1, 1), GradientG(1, 2), GradientGT(1, 3), VelocityDAtEnd(1),
        GradientG(2, 0), GradientG(2, 1), GradientG(2, 2), GradientGT(2, 3), VelocityDAtEnd(2);

    // check order
    resultMatrix /= SurfacePosition[4];
    return resultMatrix;
}

mogEigenVectorValues RecirculationSurfaceUtils::ComputeEigenVectorsValuesAt(const Vector5 SurfacePosition, const bool DoNotFilter)
{
    mogEigenVectorValues Result;
    Result.atPosition = SurfacePosition;

    Eigen::Vector3d integratedPosition = FlowData.IntegrateOverTime(SurfacePosition);

    Eigen::Matrix<double, 3, 5> SurfaceMatrix = CreateSurfaceMatrixM(SurfacePosition, integratedPosition);
    Eigen::Matrix<double, 5, 5> GrowthMatrix = SurfaceMatrix.transpose() * SurfaceMatrix;

    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> RecircEigenSolver(GrowthMatrix);
    //Eigen::EigenSolver<Eigen::MatrixXd> RecircEigenSolver(GrowthMatrix);

    auto EigenVectors = RecircEigenSolver.eigenvectors();
    auto EigenValues = RecircEigenSolver.eigenvalues();

    // find EigenVectors with Eigenvalues zero
    for (int i = 0; i < EigenValues.size(); i++) {
        auto eigenVal = EigenValues(i);
        if (DoNotFilter || fabs(eigenVal) < 1.0e-10) {
            Result.eigenVectors.push_back((Vector5() << EigenVectors.col(i)).finished());
            Result.eigenVectors.back().normalize(); // normalize each of them.
        }
    }

    Result.eigenValues = EigenValues;

    return Result;
}

const std::vector<mogDistanceResult> RecirculationSurfaceUtils::ComputeAllDistancesFor(const Vector5 StartPosition, const Eigen::Vector2d& MaxTimeTau, const unsigned int TimeSteps)
{
    // Resolution, (timesteps)
    // min max. 
    assert(StartPosition[4] != 0.0 && "Integration time is 0, need a different minimum Tau / Integration time.");

    double CellSizeT = fabs(MaxTimeTau[0] - StartPosition[3]) / static_cast<double>(TimeSteps);
    double CellSizeTau = fabs(MaxTimeTau[1] - StartPosition[4]) / static_cast<double>(TimeSteps);

    // for each t (stepsize)
    // integrate for each Tau (regular )
    // save distance 
    std::vector<mogDistanceResult> result;

    Vector5 StartVector;
    for (double tau = 0; tau < TimeSteps; tau++) {
        for (double t = 0; t < TimeSteps; t++) {
            // integrate from t over tau, 
            // compute distance from final to start position, 
            // Save / print distance
            StartVector = (Vector5() << StartPosition[0], StartPosition[1], StartPosition[2], (StartPosition[3] + t * CellSizeT), (StartPosition[4] + tau * CellSizeTau)).finished();
            
            // distance can be calculated from DistanceVector on export
            result.push_back(mogDistanceResult(StartVector, FlowData.IntegrateOverTime(StartVector) - Eigen::Vector3d(StartPosition[0], StartPosition[1], StartPosition[2])));
        }
    }

    return result;
}

void RecirculationSurfaceUtils::StartDistanceCalculation(const Vector5 MinVector, const Vector5 MaxVector, const unsigned int TimeSteps, const bool TauComparison)
{
    Eigen::Vector3d CellSizes;
    //check(GradientFieldResolution.GetMin() <= 0 && "Resolution axis invalid (needs to be greater than 0)");

    CellSizes[0] = fabs(MaxVector[0] - MinVector[0]) / (double)parameters.FieldResolution.x();
    CellSizes[1] = fabs(MaxVector[1] - MinVector[1]) / (double)parameters.FieldResolution.y();
    CellSizes[2] = fabs(MaxVector[2] - MinVector[2]) / (double)parameters.FieldResolution.z();
    auto CellSizet = fabs(MaxVector[3] - MinVector[3]) / (double)TimeSteps; // t stepsize
    if (CellSizet == 0.0)
        CellSizet = 1.0;
    // Tau stepsize varies with t.
    double endTau = MaxVector[4];

    MinimumFutureDistances.clear(); // finish everything first. // should call destructor

                                    // save each distance at the correct location
    for (int z = 0; z < parameters.FieldResolution.z(); z++) {
        for (int y = 0; y < parameters.FieldResolution.y(); y++) {
            for (int x = 0; x < parameters.FieldResolution.x(); x++) {
                
                Eigen::Vector3d StartVector3 = Eigen::Vector3d(MinVector[0] + x*CellSizes[0], MinVector[1] + y*CellSizes[1], MinVector[2] + z*CellSizes[2]);
				std::clog << "Adding task for pos. (" << StartVector3.x() << ", " << StartVector3.y() << ", " << StartVector3.z() << ") " << std::endl; // flush immediately
                MinimumFutureDistances.push_back(std::async(std::launch::async, [StartVector3, CellSizet, TauComparison, endTau, this]
                {
                    return this->ComputeMinimumDistanceWithin(StartVector3, CellSizet, endTau, TauComparison);
                })); //async on VS should use thread pools or similar (not start more threads than hardware supports / reuse threads) 
                // GCC / boost don't seem to implement it this way? (old GCC doesn't even start a new thread. :/ )
                // http://stackoverflow.com/questions/15666443/which-stdasync-implementations-use-thread-pools
                
                std::clog << "Number of Tasks: " << MinimumFutureDistances.size() << std::endl;
            }
        }
    }
}

mogDistanceResult RecirculationSurfaceUtils::ComputeMinimumDistanceWithin(const Eigen::Vector3d StartVector3, const double CellSize, double endTau, const bool TauComparison)
{
    // Find smallest distance for varying t and Tau. 
    // save distance to file
    // all the times from minimum t to maximum t
    mogDistanceResult curResultData = mogDistanceResult();

    // only need to go forward and backward and check all points in solution.
    if (endTau == 0.0) {
        endTau = fabs(parameters.MaximumBounds[3] - parameters.MinimumBounds[3]);
    }

    assert(endTau != 0.0 && "tau is zero. Is minimum and maximum starttime (w) the same?");

    //Vector5 minimumVector; // save minimum Vector as well?
    // got from minimum start time to maximum start time.
    for (int t = 0; (parameters.MinimumBounds[3] + t*CellSize) <= parameters.MaximumBounds[3]; t++) {

        Eigen::Vector4d StartVector = Eigen::Vector4d(StartVector3[0], StartVector3[1], StartVector3[2], parameters.MinimumBounds[3] + static_cast<double>(t)*CellSize);
        Solution3D currentSolution;

        //FORWARD

        // check if full integration is already minimum.
        double currentDistance;
        Eigen::Vector3d currentDistanceVector = FlowData.IntegrateOverTime(StartVector, endTau, CellSize*0.1, &currentSolution) - Eigen::Vector3d(StartVector3[0], StartVector3[1], StartVector3[2]);

        double lastDistance = currentDistanceVector.norm();		// check all distances on integrated line for optimal minimum
        curResultData.DistanceVector = currentDistanceVector;
        curResultData.MinPosition = (Vector5() << StartVector, endTau).finished(); // this will be used if no recirculation can be found.
		// TODO: Normal? -> EigenVectors?
        if (TauComparison) // we want to compare with other taus (one field only)
        {
            for (unsigned int indexSol = 1; indexSol < currentSolution.size(); ++indexSol) 
            {
                // get next distance and distance vector from solution (all points/positions of integration)
                currentDistanceVector = (currentSolution.y()[indexSol] - StartVector3);
                currentDistance = currentDistanceVector.norm();

                //first check if actually bending back to origin, then compare distance. (if the distance is only growing, it's not a recirculation.)
                if (lastDistance > currentDistance && currentDistance  < curResultData.DistanceVector.norm()) {
                    curResultData.DistanceVector = currentDistanceVector;
                    curResultData.MinPosition = (Vector5() << StartVector[0], StartVector[1], StartVector[2], StartVector[3], currentSolution.t()[indexSol - 1] - StartVector[3]).finished();
                }

                // save new distance for next iteration
                lastDistance = currentDistance;
            }
        }
        else {

            curResultData.myIntegration = currentSolution; // copy, NOTE: make more efficient if needed (unique ptr)
        }


        // currentSolution should be removed / freed here anyway.
    }

    return curResultData;
}
