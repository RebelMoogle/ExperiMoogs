#pragma once

#include "ofMain.h"
#include <ofxGui.h>
#include "TypesAndStructs.h"
#include "AnalyticField.h"
#include "RecirculationSurfaceUtils.h"
#include <future>


class ofPolyline;
class ofxTextField;

class ofApp : public ofBaseApp{

	public:
		void setup();
		void update();
		void draw();

		void keyPressed(int key);
		void keyReleased(int key);
		void mouseMoved(int x, int y );
		void mouseDragged(int x, int y, int button);
		void mousePressed(int x, int y, int button);
		void mouseReleased(int x, int y, int button);
		void mouseEntered(int x, int y);
		void mouseExited(int x, int y);
		void windowResized(int w, int h);
		void dragEvent(ofDragInfo dragInfo);
		void gotMessage(ofMessage msg);
		
private:

    ofShader IllumLines;
	// load csv eigendata button
	ofxPanel gui;
    ofxLabel CamPos;
	// button for loading exported surface data including first two eigenvectors and eigenvalues.
	ofxButton computePathLine;
	ofxButton computeDistances;
	ofxToggle RenderDistances;
	ofxButton reloadShader;
    //ofxButton loadFlowData;
    
    // menu to choose analytic field


    

	ofEasyCam cam;

    ofLight sun;
	ofSpherePrimitive testSphere;
    std::unique_ptr<AnalyticField> Flow;
	std::unique_ptr<RecirculationSurfaceUtils> FlowTools;
	std::vector<ofPolyline> Pathlines;
	std::vector<ofMesh> LineMeshes;
	std::vector < std::future< PathlineTimes> > futurePathLines;

    // computes and adds pathline
    void OnComputePathLinePress();
	void OnComputeDistances();
	void OnShaderReload();



	PathlineTimes ComputeAndAddPathline(const double x, const double y, const double z, const double t = 1.0);
	
	const std::vector<mogDistanceResult>* MinimumDistances = nullptr;



};
