#pragma once
#include <future>
#include "TypesAndStructs.h"
#include "AnalyticField.h"

struct RecirculationParameters {

    RecirculationParameters() 
        :   InitialIntegrationTime(3.0), 
            VerticesPerAxis(11), 
            InitialStepSize(0.001), 
            MaxError(0.00001), 
            FieldStartTimeStepSize(0.1), 
            ComputeSeedPoint(0), 
            FieldResolution(Eigen::Vector3i(32, 32, 32)), 
            MinimumBounds(Eigen::Vector4d(-1, -1, -1, -1)), 
            MaximumBounds(Eigen::Vector4d(1, 1, 1, 1))
    {}
    /* sort these PARAMETERS -> save in extra options struct? Or via UI?*/
    /** How long should the integration last, be computed for? (Also used for Distance Calculation as Max integration time., minimum Integration time depends on StepSize. */
    double InitialIntegrationTime;

    /** Number of data points per Axis, odd numbers are best. */
    UINT32 VerticesPerAxis;

    /** What step size should be used for the integration? */
    double InitialStepSize;

    /** What is the maximum allowable error for finding points on the surface? */
    double MaxError;

    //int FieldsWithNoTauComparison = 0;

    /* stepSize for different starttimes from minimum to maximum t (w) */
    double FieldStartTimeStepSize;

    /*	Computes number of seedpoints to use for Surface calculation. Will try to find best local minima for each point and choose the one with the least distance.
    0 means choose the given seed point and ignore minima.*/
    UINT32 ComputeSeedPoint;

    /*Resolution of the Field grid. Needs to be larger than 0*/
    Eigen::Vector3i FieldResolution;

    /*Minimum Values for the bounds for the minimum distance and Gradient field calculation (3D Space coordinates + t as 4th. dimension)*/
    Eigen::Vector4d MinimumBounds;

    /*Maximum Values for the bounds for the minimum distance and Gradient field calculation (3D Space coordinates + t as 4th. dimension)*/
    Eigen::Vector4d MaximumBounds;
};

/** Contains members and functions necessary for calculating Recirculation Surfaces */
class RecirculationSurfaceUtils {
public:
    //TODO setting parameters
    RecirculationSurfaceUtils(const AnalyticField& flowData, const RecirculationParameters& recircOptions = RecirculationParameters()) : FlowData(flowData), parameters(recircOptions)
    {};
    ~RecirculationSurfaceUtils();

    const mogDistanceResult FindSurfacePositionFrom(Vector5 StartPosition, double StepSize, double MinimumStep = 0.0, const int maximumIterations = 500);

private:
    const Vector5 ComputeGradientDescentToSurface(const Vector5& StartVector, const Eigen::Vector3d IntegratedPosition, const double StepSize = 0.001) const;

    const Eigen::Matrix<double, 3, 5> CreateSurfaceMatrixM(const Vector5 SurfacePosition, const Eigen::Vector3d IntegratedPosition, double StepSize = 0.001) const;
    mogEigenVectorValues ComputeEigenVectorsValuesAt(const Vector5 SurfacePosition, const bool DoNotFilter = false);

    /*	Calculate all distances for different start times and integration periods for a single position.
    *	StartPosition contains Position and minimum starttime and integration period,
    *	MaxTimeTau holds the maximum Starttime and Integration period #
    *   Uses Timesteps as resolution*/
    const std::vector<mogDistanceResult> ComputeAllDistancesFor(const Vector5 StartPosition, const Eigen::Vector2d MaxTimeTau);

    void StartDistanceCalculationThreaded(const Vector5 MinVector, const Vector5 MaxVector, const bool TauComparison = true);
    mogDistanceResult ComputeMinimumDistanceWithin(const VC::math::VecN<double, 3> StartVector3, const double CellSize, double endTau, const bool TauComparison = true);

    std::vector<mogSurfaceResult> SurfaceVertices;
    std::vector<std::future<mogDistanceResult>> MinimumFutureDistances; // Minimum Distances of the whole 5d field at every 3d point.
    std::vector<mogDistanceResult> MinimumDistances; //move to extra class / funktor?

    const AnalyticField& FlowData;

    RecirculationParameters parameters;

};

