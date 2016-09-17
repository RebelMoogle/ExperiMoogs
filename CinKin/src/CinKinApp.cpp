#include "cinder/app/App.h"
#include "cinder/gl/gl.h"

#include "Kinect2.h"
#include <map>

#include "LineTracker.h"


class CinKinApp : public ci::app::App
{

public:
	CinKinApp();

	void update() override;
	void draw() override;
private:
	ci::Channel16uRef mChannel;
	Kinect2::DeviceRef mDevice;
	std::vector<Kinect2::Body> mBodies;
	ci::geom::Sphere mSphere = ci::geom::Sphere(ci::Sphere(glm::vec3(0., 0., 0.), 10.0));
	glm::vec3 position;
	ci::CameraPersp			mCam;
	std::map<uint64_t, std::vector<LineTracker>> allLines;

	std::vector < JointType> trackJoints;

};

#include "cinder/app/RendererGl.h"

CinKinApp::CinKinApp()
{	
	mDevice = Kinect2::Device::create();
	mDevice->start();
	mDevice->enableHandTracking();
	mDevice->connectDepthEventHandler( [ & ]( const Kinect2::DepthFrame& frame )
	{
		mChannel = frame.getChannel();
	} );
	mDevice->connectBodyEventHandler([&](const Kinect2::BodyFrame& frame) {
	
		mBodies = frame.getBodies();
	
	});

	mCam.lookAt(glm::vec3(0, 0, -2), glm::vec3(0));

	ci::gl::enableDepthWrite();
	ci::gl::enableDepthRead();

	trackJoints.push_back(JointType_HandTipRight);
	trackJoints.push_back(JointType_HandTipLeft);
	trackJoints.push_back(JointType_FootLeft);
	trackJoints.push_back(JointType_FootRight);
	trackJoints.push_back(JointType_Head);
}

void CinKinApp::update()
{
	// Update line trackers
	for (std::vector<Kinect2::Body>::iterator itBody = mBodies.begin(); itBody != mBodies.end(); ++itBody)
	{
		if (itBody->isTracked() || itBody->isEngaged())
		{
			uint64_t bodyID = itBody->getId();
			std::map<uint64_t, std::vector<LineTracker>>::iterator itLine = allLines.find(bodyID);
			if (itLine == allLines.end())
			{
				// linetrackers not found, create
				std::pair<std::map<uint64_t, std::vector<LineTracker>>::iterator, bool> itPair = allLines.insert({ bodyID, std::vector<LineTracker>() });

				assert(itPair.second);
				itLine = itPair.first;

			}

			std::vector<LineTracker>& bodyLines = itLine->second;
			std::map<JointType, Kinect2::Body::Joint>  zeJoints = itBody->getJointMap();


			while (bodyLines.size() < 5)
			{
				bodyLines.push_back(LineTracker());
			}

			for (LineTracker curTracker : bodyLines)
			{
				curTracker.update();
			}

			int JointCount = 0;
			for( JointType trackType : trackJoints)
			{

				bodyLines[JointCount].addPoint(zeJoints[trackType].getPosition());
				++JointCount;
			}
		}
	}
}

void CinKinApp::draw()
{
	/*if ( mChannel ) {
		const ci::gl::TextureRef tex = ci::gl::Texture::create( *Kinect2::channel16To8( mChannel ) );
		ci::gl::draw( tex, tex->getBounds(), getWindowBounds() );
	}*/

	ci::gl::clear(cinder::ColorAT<float>(0, 0, 0));

	ci::gl::setMatrices(mCam);

	for each (std::pair<uint64_t, std::vector<LineTracker>> bodyLines in allLines)
	{
		for each (LineTracker currentLine in bodyLines.second)
		{
			currentLine.draw();
		}
	}


	


	
	//ci::gl::drawColorCube(glm::vec3(0, 0, 0), glm::vec3(1,1,1));
	
}

CINDER_APP( CinKinApp, ci::app::RendererGl )
	
