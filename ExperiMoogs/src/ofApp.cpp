#include "ofApp.h"

//--------------------------------------------------------------
void ofApp::setup(){
	ofSetVerticalSync(true);

	RandomColor.addListener(this, &ofApp::ChangeToRandomColor);

	gui.setup();
	gui.add(RandomColor.setup("Random Color"));
	gui.add(ChangeTime.setup("color changing period", 1, 0.1, 10));
	gui.add(ChooseColor.setup("color", ofColor(100, 100, 140), ofColor(0, 0, 0), ofColor(255, 255, 255)));



}

//--------------------------------------------------------------
void ofApp::update(){

	auto delta = ofGetLastFrameTime();

	if(remainingTime < ChangeTime)
	{
		currentColor.lerp(nextColor, delta / ChangeTime);
		remainingTime += delta; // should be in seconds?
	}
	else
	{
		ChangeToRandomColor();
	}

}

//--------------------------------------------------------------
void ofApp::draw(){
	ofBackgroundGradient(currentColor, ChooseColor);

	gui.draw();
}

//--------------------------------------------------------------
void ofApp::keyPressed(int key){

}

//--------------------------------------------------------------
void ofApp::keyReleased(int key){

}

//--------------------------------------------------------------
void ofApp::mouseMoved(int x, int y ){

}

//--------------------------------------------------------------
void ofApp::mouseDragged(int x, int y, int button){

}

//--------------------------------------------------------------
void ofApp::mousePressed(int x, int y, int button){

}

//--------------------------------------------------------------
void ofApp::mouseReleased(int x, int y, int button){

}

//--------------------------------------------------------------
void ofApp::mouseEntered(int x, int y){

}

//--------------------------------------------------------------
void ofApp::mouseExited(int x, int y){

}

//--------------------------------------------------------------
void ofApp::windowResized(int w, int h){

}

//--------------------------------------------------------------
void ofApp::gotMessage(ofMessage msg){

}

void ofApp::ChangeToRandomColor()
{
	nextColor.set(ofRandom(0, 255), ofRandom(0, 255), ofRandom(0, 255));
	remainingTime = 0;
}

//--------------------------------------------------------------
void ofApp::dragEvent(ofDragInfo dragInfo){ 

}
