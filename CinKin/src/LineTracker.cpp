#include "LineTracker.h"

LineTracker::~LineTracker()
{
}

void LineTracker::addPoint(const glm::vec3 & position)
{
		trackedLine.push_back(position);

		int tooMany = trackedLine.size() - m_MaxPoints;
		if (tooMany > 0)
		{
			std::vector<glm::vec3>& linePoints = trackedLine.getPoints();
			linePoints.erase(linePoints.begin(), linePoints.begin() + tooMany);
		}

}

void LineTracker::draw()
{
	ci::gl::color(lineColor);
	glm::vec3 position = trackedLine.getPoints().back();
	//ci::gl::drawSphere(position, 0.1);
	ci::gl::draw(trackedLine);
	//ci::gl::drawElements()
	// need shader for color? Vertex Buffer probably...  
}

void LineTracker::update()
{
	std::vector<glm::vec3>& linePoints = trackedLine.getPoints();
	if(!linePoints.empty() )
	{
		linePoints.erase(linePoints.begin());
	}
}
