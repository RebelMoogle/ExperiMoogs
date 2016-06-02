#pragma once
#include <vector>
#include "ofMain.h"

//reserves memory for number of datapoints along a unitcircle. 
// get data by angle.
template<typename DataType>
static bool setUnitCircleData(std::vector<DataType>& unitCircleData, double angle, DataType dataToSet)
{
	double correctedAngle = angle;
	// make sure angle is below 2pi
	if(fabs(correctedAngle) > M_TWO_PI)
	{
		correctedAngle = modf(fabs(correctedAngle), M_TWO_PI);
	}

	// make sure negative angles get handled directly.
	if(angle < 0.0)
	{
		correctedAngle = M_TWO_PI - correctedAngle;
	}
	// TODO: any exceptions go here
	
	// calculate index in array
	double stepSize = M_TWO_PI / static_cast<double>(unitCircleData.size); // compiler should optimize this away.
	unsigned index = correctedAngle / stepSize; // Angle could be negative, absolute value. TODO: Lerp, instead of convert to int. 
	
	unitCircleData[index] = dataToSet;
	// TODO set with gradient to sourrounding datapoints as well.
	while(stepSize * fallOff > )
	// (stepSize * falloffperStep)
	
	return true;
}

// note: this could also work as static function: set( vector&, angle, Data); get(vector&) // number of DataPoints is implicitly given by length of vector.



// Markov Chain Vector Graphics. 
// User input graphics. - record likely hood of vector following another vector. 
// recorded vector data: angle to up vector ( zenith ), (length, distance from origin of graphic)
// save likelyhood. (start with angle first, to test. use unit sphere angle? sphherical angles?
// info for new output vector angle to up vector, length
//

// record via events? (or look up observer pattern)
