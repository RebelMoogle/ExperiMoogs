#include "ofApp.h" 

//-------------------------------------------------------------- 
void ofApp::setup()
{


    gui.setup();

    // TODO: setup gui 
    // add buttons for starting and displaying calculation 
    // add parameter settings -> enable change of parameters! 
    // TODO: might have to make FlowData changeable. (pure pointer, have instance in main app).  

    cam.setAutoDistance(false);
    cam.enableMouseInput();
    cam.setNearClip(0.000001);
    cam.setFarClip(1000.0);
    cam.setPosition(ofVec3f(0, 0, -1));

    glEnable(GL_POINT_SMOOTH); // use circular points instead of square points 
    glPointSize(3); // make the points bigger 
}

//-------------------------------------------------------------- 
void ofApp::update()
{

}

//-------------------------------------------------------------- 
void ofApp::draw()
{
    ofBackgroundGradient(ofColor(64), ofColor(0));

    ofSetColor(255);
    cam.begin();

    ofSetColor(ofColor::white);
    // TODO: render surface 
    // TODO: display distances 

    cam.end();


    ofSetColor(ofColor::gray);

    ofNoFill();
    ofSetColor(ofColor::yellow);
    ofSetLineWidth(2);
    ofDrawCircle(cam.worldToScreen(selectedVert), 4);
    ofSetLineWidth(1);

    ofVec2f offset(10, -10);

    // TODO: information about highlighted vertex 


    gui.draw();
}

//-------------------------------------------------------------- 
void ofApp::keyPressed(int key)
{

}

//-------------------------------------------------------------- 
void ofApp::keyReleased(int key)
{

}

//-------------------------------------------------------------- 
void ofApp::mouseMoved(int x, int y)
{

}

//-------------------------------------------------------------- 
void ofApp::mouseDragged(int x, int y, int button)
{

}

//-------------------------------------------------------------- 
void ofApp::mousePressed(int x, int y, int button)
{


    if (button == OF_MOUSE_BUTTON_3) {
        // TODO: select vertex 
    }

}

//-------------------------------------------------------------- 
void ofApp::mouseReleased(int x, int y, int button)
{

}

//-------------------------------------------------------------- 
void ofApp::mouseEntered(int x, int y)
{

}

//-------------------------------------------------------------- 
void ofApp::mouseExited(int x, int y)
{

}

//-------------------------------------------------------------- 
void ofApp::windowResized(int w, int h)
{

}

//-------------------------------------------------------------- 
void ofApp::gotMessage(ofMessage msg)
{

}

//-------------------------------------------------------------- 
void ofApp::dragEvent(ofDragInfo dragInfo)
{

}
\ No newline at end of file
Moving to boost odeint

Integrated boost::numeric::odeint
->next testing