#include "ofApp.h" 
#include "ofPolyline.h"
#include "FlowData.h"

//-------------------------------------------------------------- 
void ofApp::setup()
{
    ofSetSmoothLighting(true);

    sun.setup();
    sun.setDirectional();
    sun.setAmbientColor(ofFloatColor(0.1, 0.1, 0.1));
    sun.setDiffuseColor(ofFloatColor(1, 1, 1));
    sun.setSpecularColor(ofFloatColor(1, 1, 1));
    sun.setOrientation(ofVec3f(15, 30, 0));

    // shininess is a value between 0 - 128, 128 being the most shiny //
    material.setColors(ofFloatColor::yellow, ofFloatColor::lightYellow, ofFloatColor::white, ofFloatColor::black);
    material.setShininess(120);
    // the light highlight of the material //

    // AnalyticField(std::string fieldName, std::function<Eigen::Vector3d(const Eigen::Vector3d& x, Eigen::Vector3d& dxdt, const double t)> analyticFormula);
    MyAnalyticField = make_unique<AnalyticField>("DoubleGyre3D", FlowData::DoubleGyre3D);
    computePathLine.addListener(this, &ofApp::OnComputePathLinePress);



    gui.setup();

    // button for calculating pathline
    gui.add(computePathLine.setup("Compute Pathline"));
    gui.add(CamPos.setup("Camera Position", "Camera Position"));

    // TODO: might have to make FlowData changeable. (pure pointer, have instance in main app).  
    // choose / load FlowData button (choose from menu.) -> Later load from file. (Lua or Julia? )
    

    // draw/calculate pathline button

    cam.setAutoDistance(false);
    cam.enableMouseInput();
    cam.setNearClip(0.000001);
    cam.setFarClip(100.0);
    cam.setPosition(ofVec3f(1.0, 0.3, 1.0));
    cam.lookAt(ofVec3f(0., 0., 0.), ofVec3f(0., 1., 0.));

    glEnable(GL_POINT_SMOOTH); // use circular points instead of square points 
    glPointSize(3); // make the points bigger 
}

//-------------------------------------------------------------- 
void ofApp::update()
{
    CamPos = std::to_string(cam.getPosition().x) + ", " + std::to_string(cam.getPosition().y) + ", " + std::to_string(cam.getPosition().z);

    // add completed path line calculations
    if (!futurePathLines.empty())
    {
        if (futurePathLines.front().valid()) {
            PathLines.push_back(futurePathLines.front().get());
        }
        else {
            // remove if invalid.
            futurePathLines.erase(futurePathLines.begin());
        }
    }


}

//-------------------------------------------------------------- 
void ofApp::draw()
{
    ofBackgroundGradient(ofColor(64), ofColor(0));

    ofEnableLighting();
    sun.enable();
    cam.begin();

    ofSetColor(ofColor::white);
    // TODO: render surface 
    // TODO: display distances 

    ofDrawGrid(0.1f, 10, true);

    ofSetColor(ofColor::yellow);
    ofSetLineWidth(3);
    //material.begin();

    for (auto curLine : PathLines) {
        curLine.draw();
    }
    //material.end();
    ofSetLineWidth(1);

    cam.end();
    sun.disable();
    ofDisableLighting();

    ofSetColor(ofColor::white);


    

    //ofVec2f offset(10, -10);

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

void ofApp::OnComputePathLinePress()
{
    // TODO parameter?
    const double stepSize = 0.1;

    // create 10 * 10 * 10 pathlines, with 0.1 stepsize

    for (int x = 0; x < 10; ++x) {
        for (int y = 0; y < 10; ++y) {
            for (int z = 0; z < 10; ++z) {
                futurePathLines.push_back(std::async(std::launch::async, [x, y, z, stepSize, this]
                {
                    return this->ComputeAndAddPathline(x * stepSize, y * stepSize, z * stepSize);
                }));
            } // for z
        } // for y
    } // for x
}

ofPolyline ofApp::ComputeAndAddPathline(const double x, const double y, const double z, const double t)
{
    if (!MyAnalyticField) {
        throw std::exception("Analytic Flow Data not set!");
    }
    ofPolyline pathline;
    MyAnalyticField->ComputePathLineAt(Eigen::Vector4d(x, y, z, t), pathline);
    return pathline;

}
             

//-------------------------------------------------------------- 
void ofApp::dragEvent(ofDragInfo dragInfo)
{

}