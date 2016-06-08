#pragma once

#include "ofMain.h"
#include <ofxGui.h>
#include "TypesAndStructs.h"
#include "AnalyticField.h"
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
		
	// load csv eigendata button
	ofxPanel gui;
    ofxLabel CamPos;
	// button for loading exported surface data including first two eigenvectors and eigenvalues.
	ofxButton computePathLine;
    //ofxButton loadFlowData;
    
    // menu to choose analytic field


    

	ofEasyCam cam;

    ofLight sun;
    ofMaterial material;

	private:
        std::unique_ptr<AnalyticField> MyAnalyticField;
        std::vector<ofPolyline> PathLines;
        std::vector<std::future<ofPolyline>> futurePathLines;

        // computes and adds pathline
        void OnComputePathLinePress();

        ofPolyline ComputeAndAddPathline(const double x, const double y, const double z, const double t = 1.0);




};
