# ----

Prerequesites. 
*OpenFrameworks
Currently using OpenFrameworks, (http://openframeworks.cc/)
Will need recompile and linking for your platform. 
Everything else is header only. 

*Boost::numeric::odeint 
odeint required as well. (Either install boost completely, or get just odeint sources: https://headmyshoulder.github.io/odeint-v2/)
odeint is header-only, but seems to require some lib files anyway. :/

*Eigen
Eigen is also header only and included in the repo (under ThirdParty)


Settings. 
Include Directories: 
*(PROJECT_DIR)\ThirdParty
*BOOST ROOT DIRECTORY (BOOST_ROOT path variable?

Library Directories: (Linker)
*(BOOST_ROOT)\<PLATFORM LIBRARY FILES> --> might need to download these from the boost site or compile yourself.

(vclibs currently not in use, since it requires compilation and is not trivial to set up and get working (on Win). :/ )

Setup: 

Easiest way is the following: 
1. clone this repo
2. Download Openframeworks for your platform.
3. Use Projectgenerator to create project files. (currently only using OFxGui addon)
4. Add include directories for Eigen and Boost in IDE. 

You can of course add openframeworks by hand. 

At some point I might create an automated version, either creating ofxaddons or using cmake or similar. 
Might also remove openframeworks dependency and replace with imgui and pure opengl, but that will require a bit of extra work. 
(mostly reimplementing functionality in OpenGL)