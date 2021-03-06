#pragma once
#include <future>
#include "AnalyticField.h"
#include "TypesAndStructs.h"

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
    unsigned int VerticesPerAxis;

    /** What step size should be used for the integration? */
    double InitialStepSize;

    /** What is the maximum allowable error for finding points on the surface? */
    double MaxError;

    //int FieldsWithNoTauComparison = 0;

    /* stepSize for different starttimes from minimum to maximum t (w) */
    double FieldStartTimeStepSize;

    /*	Computes number of seedpoints to use for Surface calculation. Will try to find best local minima for each point and choose the one with the least distance.
    0 means choose the given seed point and ignore minima.*/
    unsigned int ComputeSeedPoint;

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
    RecirculationSurfaceUtils(const AnalyticField& flowData, const RecirculationParameters& recircOptions = RecirculationParameters()) : FlowData(flowData), parameters(recircOptions)
    {};
    ~RecirculationSurfaceUtils();

    const mogDistanceResult FindSurfacePositionFrom(Vector5 StartPosition, double StepSize, double MinimumStep = 0.0, const int maximumIterations = 500);

    /* gets all calculated distances. 
    * if the calculation is not done yet, it will clear the old saved results and 
    fill MinimumDistances with the new calculation results. (waiting for them to finish)
    * throws an exception if calculations have not been started.**/
    const std::vector<mogDistanceResult>& GetDistances();

    /* Start Distance Calculation in threads
    * uses std::future for parralel computing. 
    * gcc: needs a current version of gcc (5+) to run correctly.  **/
    void StartDistanceCalculation(const ofVec4f MinVector, const ofVec4f MaxVector, const float tau, const unsigned int TimeSteps = 11, const bool TauComparison = true);

private:

    const Vector5 ComputeGradientDescentToSurface(const Vector5 StartVector, const Eigen::Vector3d IntegratedPosition, const double StepSize = 0.001) const;

    const Eigen::Matrix<double, 3, 5> CreateSurfaceMatrixM(const Vector5 SurfacePosition, const Eigen::Vector3d IntegratedPosition, double StepSize = 0.001) const;
    mogEigenVectorValues ComputeEigenVectorsValuesAt(const Vector5 SurfacePosition, const bool DoNotFilter = false);

    /*	Calculate all distances for different start times and integration periods for a single position.
    *	StartPosition contains Position and minimum starttime and integration period,
    *	MaxTimeTau holds the maximum Starttime and Integration period #
    *   Uses Timesteps as resolution*/
    const std::vector<mogDistanceResult> ComputeAllDistancesFor(const Vector5 StartPosition, const Eigen::Vector2d& MaxTimeTau, const unsigned int TimeSteps);
    mogDistanceResult ComputeMinimumDistanceWithin(const Eigen::Vector3d StartVector3, const double CellSize, double endTau, const bool TauComparison = true);

    std::vector<mogSurfaceResult> SurfaceVertices;
    std::vector<std::future<mogDistanceResult>> MinimumFutureDistances; // Minimum Distances of the whole 5d field at every 3d point.
    std::vector<mogDistanceResult> MinimumDistances; //move to extra class / funktor?

    const AnalyticField& FlowData;

    RecirculationParameters parameters;
    
};

