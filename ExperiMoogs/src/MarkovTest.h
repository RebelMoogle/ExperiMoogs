#pragma once
#include <map>
#include "ofMain.h"

class MarkovTest
{
public:
	MarkovTest();
	~MarkovTest();

private:
	// save probabilities. increase probability for each vector. 
	// need to save probability of all angles per vector. 
	// either a probability matrix, or a map in a map, all angle-probabilities per angle.
	// matrix would constantly change as new angles get added. when new angle gets added, does it need to get added to all existing? NO-> just the observed one. 
	// Problem: when recording a new vector angle I need to search by angle, but when accessing for application I need to search by probability. 

	// set resolution: create array/structure with enough values, between 0 - 2pi, add probablities not exactly but with a gradient to multiple values.  (falloff range)
	// vector of probability-data pairs. with numValues
	// need to save all angles per angle, not just one. 
	std::vector<std::vector<double>> angleProbabilityPairs; // create new unitcircle data type. saves arbitrary data for each data point at a resolution and handles getter and setter methods. 

};



// Markov Chain Vector Graphics. 
// User input graphics. - record likely hood of vector following another vector. 
// recorded vector data: angle to up vector ( zenith ), (length, distance from origin of graphic)
// save likelyhood. (start with angle first, to test. use unit sphere angle? sphherical angles?
// info for new output vector angle to up vector, length
//

// record via events? (or look up observer pattern)
