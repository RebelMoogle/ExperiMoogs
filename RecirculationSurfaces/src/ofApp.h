#pragma once

#include "ofMain.h"
#include <ofxGui.h>
#include "TypesAndStructs.h"
#include "AnalyticField.h"


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
	// button for loading exported surface data including first two eigenvectors and eigenvalues.
	ofxButton computePathLine;
    //ofxButton loadFlowData;
    
    // menu to choose analytic field


    

	ofEasyCam cam;

	private:
        std::unique_ptr<AnalyticField> MyAnalyticField;
        std::vector<ofPolyline> PathLines;
        
        // computes and adds pathline
        void ComputeAndAddPathline();
};
