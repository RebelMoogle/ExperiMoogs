#pragma once

#include "ofMain.h"
#include "ofxGui.h"

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

		void ChangeToRandomColor();
		
		ofxPanel gui;
		ofxFloatSlider ChangeTime;
		ofxButton RandomColor;
		ofxColorSlider ChooseColor;

	private:

		ofColor currentColor = ofColor::black;
		ofColor nextColor = ofColor::black;
		float remainingTime = 0;
};
