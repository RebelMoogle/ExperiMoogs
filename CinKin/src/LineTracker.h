#pragma once

#include "cinder\Cinder.h"
#include "cinder\Rand.h"
#include "cinder\Color.h"
#include "cinder\PolyLine.h"
#include "cinder\gl\gl.h"

class LineTracker
{
public:
	LineTracker()
	{
		lineColor = ci::Color(ci::randFloat(), ci::randFloat(), ci::randFloat());
	};
	LineTracker(unsigned int maxPoints) : m_MaxPoints(maxPoints) 
	{
		lineColor = ci::Color(ci::randFloat(), ci::randFloat(), ci::randFloat());
	};
	~LineTracker();

	const unsigned int m_MaxPoints = 100;

	void addPoint(const glm::vec3& position);

	void addPointColor(const glm::vec3& position, const ci::Color& color)
	{
		trackedLine.push_back(position);
	}

	void draw();
	void update();


private: 
	ci::PolyLine3 trackedLine;
	ci::Color lineColor;
};

