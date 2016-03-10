#include "ofApp.h"

//--------------------------------------------------------------
void ofApp::setup(){
	
	LoadEigenCSV.addListener(this, &ofApp::LoadSurfaceEigenData);
	
	gui.setup();
	gui.add(LoadEigenCSV.setup("Load Surface / EigenVector CSV Data")); 
	gui.add(DistanceThreshold.setup("Distance threshold", 0.1, 0.0, 1));
	gui.add(ArrowScale.setup("Arrow scale", 0.01, 0.0, 0.5));


	cam.setAutoDistance(false);
	cam.enableMouseInput();
	cam.setNearClip(0.000001);
	cam.setFarClip(1000.0);
	cam.setPosition(ofVec3f(0, 0, -1));

	glEnable(GL_POINT_SMOOTH); // use circular points instead of square points
	glPointSize(3); // make the points bigger
}

//--------------------------------------------------------------
void ofApp::update(){

}

//--------------------------------------------------------------
void ofApp::draw(){
	ofBackgroundGradient(ofColor(64), ofColor(0));

	ofSetColor(255);
	cam.begin();

	ofSetColor(ofColor::white);
	if (surfaceLoaded)
	{
		SurfaceData.surfaceMesh.drawVertices();

		for (int index : nearestIndices)
		{
			ofVec3f vertex = SurfaceData.surfaceMesh.getVertex(index);
			Vector5 Eigen1 = SurfaceData.FirstZeroEigenVector[index];
			Vector5 Eigen2 = SurfaceData.SecondZeroEigenVector[index];

			ofSetColor(ofColor::purple);
			ofDrawArrow(vertex,  vertex + ofVec3f(Eigen1.x(), Eigen1.y(), Eigen1.z()).normalize() * ArrowScale, ArrowScale*0.1);

			ofSetColor(ofColor::orange);
			ofDrawArrow(vertex, vertex + ofVec3f(Eigen2.x(), Eigen2.y(), Eigen2.z()).normalize() * ArrowScale, ArrowScale*0.1);
		}
	}

	cam.end();

	
	ofSetColor(ofColor::gray);

	ofNoFill();
	ofSetColor(ofColor::yellow);
	ofSetLineWidth(2);
	ofDrawCircle(cam.worldToScreen(selectedVert), 4);
	ofSetLineWidth(1);

	ofVec2f offset(10, -10);


	ofDrawBitmapStringHighlight(ofToString(selectedVert) + "\n ", ofVec2f(mouseX, mouseY) + offset);


	// todo print eigenvalues of selected point
	// highlight selected point. // see point picker example

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
	
	
	if(button == OF_MOUSE_BUTTON_3)
	{
		if (surfaceLoaded)
		{
			int n = SurfaceData.surfaceMesh.getNumVertices();
			float nearestDistance = 0;
			ofVec2f mouse(x, y);
			nearestIndices.clear();

			// find selected vertex // TODO: Improve Picking: replace with raytrace - sphere. (limit vertices to select)
			for (int i = 0; i < n; ++i)
			{
				ofVec3f cur = SurfaceData.surfaceMesh.getVertex(i);
				float distance = cam.worldToScreen(cur).distance(mouse);
				if (i == 0 || distance < nearestDistance) {
					nearestDistance = distance;
					selectedVert = cur;
				}
			}

			// select nearest vertices in 3D. 
			for (int i = 0; i < n; ++i)
			{
				ofVec3f cur = SurfaceData.surfaceMesh.getVertex(i);
				float distance = cur.distance(selectedVert);
				if (distance < DistanceThreshold)
				{
					nearestIndices.push_back(i);
				}
			}
			
			cam.setTarget(selectedVert);
		}
	}

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

//--------------------------------------------------------------
void ofApp::dragEvent(ofDragInfo dragInfo){ 

}


void ofApp::LoadSurfaceEigenData()
{
	// open file dialog (to find path and file)
	ofFileDialogResult openFileResult = ofSystemLoadDialog("Select Surface Eigen Data (CSV)");
	if (openFileResult.bSuccess && ofIsStringInString(openFileResult.fileName, ".csv"))
	{
		ofLogVerbose() << "File selected: " << openFileResult.fileName;
	}
	else if (!ofIsStringInString(openFileResult.fileName, ".csv"))
	{
		ofSystemAlertDialog("Incorrect file ending, please select a .csv file.");
		return;
	}
	else
	{
		ofSystemAlertDialog("Loading file cancelled or failed.");
		return;
	}

	try
	{

		SurfaceData = ImportPly::ImportEigenCloud(openFileResult.filePath);
	}
	catch (const std::runtime_error& e)
	{

		ofSystemAlertDialog(std::string("Error loading file: ") + e.what());
		return;
	}

	ofSystemAlertDialog(std::string("File loaded successfully, number of vertices: ") + ofToString(SurfaceData.surfaceMesh.getNumVertices()));
	surfaceLoaded = true;
	cam.setTarget(SurfaceData.surfaceMesh.getCentroid());
	
}