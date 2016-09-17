#include "ofApp.h"

//--------------------------------------------------------------
void ofApp::setup(){
	ofSetVerticalSync(true);
	ofSetSmoothLighting(true);
	glDisable(GL_CULL_FACE);
	ofEnableLighting();


	loadMeshButton.addListener(this, &ofApp::loadMesh);
	
	gui.setup();
	gui.add(loadMeshButton.setup("Load Mesh"));
	gui.add(colorSlider.setup("color", ofColor::aliceBlue, ofColor::black, ofColor::white));

	cam.setAutoDistance(true);
	cam.enableMouseInput();
	cam.setNearClip(0.001);
	cam.setFarClip(10000.0);
	cam.setPosition(ofVec3f(0, 0, -1));

	
	sun.setDirectional();
	sun.setAmbientColor(ofFloatColor::grey);
	sun.setAttenuation();
	sun.setDiffuseColor(ofFloatColor::white);
	sun.setPosition(10, 10, 10);
	sun.setSpecularColor(ofFloatColor::cyan);
	sun.setOrientation(ofVec3f(30, 30, 0));
	sun.setup();

}

//--------------------------------------------------------------
void ofApp::update(){

	auto delta = ofGetLastFrameTime();

}

//--------------------------------------------------------------
void ofApp::draw(){
	ofBackgroundGradient(ofColor::skyBlue, colorSlider);


	ofEnableDepthTest();
	sun.enable();
	cam.begin();
	if(loadedMesh.hasMeshes())
	{
		loadedMesh.enableMaterials();
		loadedMesh.enableNormals();
		loadedMesh.enableTextures();
		loadedMesh.enableColors();
		loadedMesh.drawFaces();
		loadedMesh.disableMaterials();
		loadedMesh.disableNormals();
		loadedMesh.disableTextures();
		loadedMesh.disableColors();
	}
	cam.end();
	sun.disable();

	ofDisableDepthTest();
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

void ofApp::loadMesh()
{
	ofFileDialogResult openFileResult = ofSystemLoadDialog("Select 3D Model to load");
	if(openFileResult.bSuccess)
	{
		ofLogVerbose() << "File Selected: " << openFileResult.fileName;
	}
	else
	{
		ofSystemAlertDialog("Loading file cancelled or failed.");
		return;
	}


	if(loadedMesh.loadModel(openFileResult.filePath))
	{
		ofLogVerbose() << "Model loaded successfully: " << loadedMesh.getMeshNames()[0];
	}
	cam.setTarget(loadedMesh.getPosition());
}

//--------------------------------------------------------------
void ofApp::dragEvent(ofDragInfo dragInfo){ 

}
