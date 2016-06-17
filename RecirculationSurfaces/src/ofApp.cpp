#include "ofApp.h" 
#include "ofPolyline.h"
#include "FlowData.h"

//-------------------------------------------------------------- 
void ofApp::setup()
{
    ofSetSmoothLighting(true);


	OnShaderReload();

	testSphere.setRadius(1.0);

    sun.setup();
    sun.setDirectional();
    sun.setAmbientColor(ofFloatColor(0.1, 0.1, 0.1));
    sun.setDiffuseColor(ofFloatColor(1, 1, 1));
    sun.setSpecularColor(ofFloatColor(1, 1, 1));
    sun.setOrientation(ofVec3f(15, 30, 0));

    // AnalyticField(std::string fieldName, std::function<Eigen::Vector3d(const Eigen::Vector3d& x, Eigen::Vector3d& dxdt, const double t)> analyticFormula
    Flow = make_unique<AnalyticField>("DoubleGyre3D", FlowData::DoubleGyre3D);
	FlowTools = std::make_unique<RecirculationSurfaceUtils>(*Flow);
    computePathLine.addListener(this, &ofApp::OnComputePathLinePress);
	reloadShader.addListener(this, &ofApp::OnShaderReload);
	computeDistances.addListener(this, &ofApp::OnComputeDistances);



    gui.setup();

    // button for calculating pathline
	gui.add(reloadShader.setup("Reload Shader"));
    gui.add(computePathLine.setup("Compute Pathline"));
	gui.add(computeDistances.setup("Compute Distances"));
	gui.add(RenderDistances.setup("Render Distances", false));
    gui.add(CamPos.setup("Camera Position", "Camera Position"));

    // TODO: might have to make FlowData changeable. (pure pointer, have instance in main App).  
    // choose / load FlowData button (choose from menu.) -> Later load from file. (Lua or Julia? )
    

    // draw/calculate Pathline button

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
			PathlineTimes pathline = futurePathLines.front().get();
			Pathlines.push_back(std::get<0>(pathline));
			std::vector<double>& times = std::get<1>(pathline);
			// convert to vertex and index buffer?
			// need vertex, normal, index?
			ofMesh tempMeshline;
			tempMeshline.disableIndices();
			tempMeshline.setMode(ofPrimitiveMode::OF_PRIMITIVE_LINE_STRIP); // draw as line strip
			tempMeshline.addVertices(Pathlines.back().getVertices());
			for (int pointIndex = 1; pointIndex < times.size(); ++pointIndex) {

				//TODO: FIX!
				//ofVec3f pointPos = GetVec3fFrom(polylineSolution.y()[pointIndex]);
				
				tempMeshline.addNormal(Pathlines.back().getNormalAtIndexInterpolated(pointIndex));
				ofVec2f ofTexCoord = ofVec2f(	static_cast<float>(pointIndex) / static_cast<float>(times.size()), 
												times[pointIndex] / times.back());
				tempMeshline.addTexCoord(ofTexCoord);
			}


			LineMeshes.push_back(tempMeshline);
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
	
	// TODO: render surface 
    // TODO: display distances 
	if (RenderDistances) 		{
		try {
			if (!MinimumDistances) 			{ MinimumDistances = &FlowTools->GetDistances();}
			// TODO: render distances.

			// TODO: either render as colored point cloud
			// or implement volume rendering. 
			// both will require conversion of data.
		}
		catch (std::exception& e) {
			std::string message = e.what();
			message += " - Distances not available, did you start the calculation? \n";
			ofLogWarning() << message;
			ofSystemAlertDialog(message);
			RenderDistances = false;
		}
	}


	ofEnableDepthTest();

    ofDrawGrid(0.1f, 10, true);

    ofSetColor(ofColor::yellow);
    ofSetLineWidth(3);

	
	IllumLines.begin();
	IllumLines.setUniform4fv("shaderColor", ofFloatColor::lightSteelBlue.v);
	IllumLines.setUniform3f("camDirection", cam.getLookAtDir());
	IllumLines.setUniform3f("camPosition", cam.getPosition());
		for (auto curLine : LineMeshes) {
			curLine.draw();
		}
	IllumLines.end();
    ofSetLineWidth(1);

	
    cam.end();
    sun.disable();
    ofDisableLighting();
	ofDisableDepthTest();

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

    for (int x = 0; x < 20; x+=2) {
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

void ofApp::OnComputeDistances()
{
	ofLogNotice() << "Starting Distance Calculation. \n";
	FlowTools->StartDistanceCalculation((Vector5()<< 0., 0., 0., 0., 0.1).finished(), (Vector5() << 2., 1., 1., 1., 5.).finished());
}

void ofApp::OnShaderReload()
{
	ofLogNotice() << "Reloading Shader.";
	IllumLines.load("simple");
	GLenum glError = glGetError();
	if (glError != GL_NO_ERROR) {
		ofLogNotice() << "Error while loading shader: " << glewGetErrorString(glError);
	}
	
}

PathlineTimes ofApp::ComputeAndAddPathline(const double x, const double y, const double z, const double t)
{
    if (!Flow) {
        throw std::exception("Analytic Flow Data not set!");
    }
	ofPolyline pathline;
	std::vector<double> times;
    Flow->ComputePathLineAt(Eigen::Vector4d(x, y, z, t), pathline, &times);
    return std::make_tuple(pathline, times);
}
             

//-------------------------------------------------------------- 
void ofApp::dragEvent(ofDragInfo dragInfo)
{

}