# ----

Prerequesites. 
*OpenFrameworks
Currently using OpenFrameworks, (http://openframeworks.cc/)
Will need recompile and linking for your platform. 
**Required OpenFrameworks Addons:
**ofxGui (included)

*Boost::numeric::odeint 
odeint required as well. (Either install boost completely, or get just odeint sources: https://headmyshoulder.github.io/odeint-v2/)
odeint is header-only, but seems to require some lib files anyway. :/

*Eigen
Eigen is also header only and included in the repo (under ThirdParty)


Settings. 
Include Directories: 
*(PROJECT_DIR)\ThirdParty
*(BOOST_ROOT) DIRECTORY (BOOST_ROOT system variable)

Library Directories: (Linker)
*(BOOST_ROOT)\<PLATFORM LIBRARY FILES> --> might need to download these from the boost site or compile yourself.

(vclibs currently not in use, since it requires compilation and is not trivial to set up and get working (on Win). :/ )

Setup: 

Easiest way is the following: 
1. clone this repo
2. Download Openframeworks for your platform.
3. Use Projectgenerator to create project files. 
	3.1. Import project folder
	3.2. add ofxGui addon
	3.3. Generate / Update project files.
4. Add include and library directories for Eigen and Boost in IDE / make manually

You can of course add openframeworks by hand. 


TODO: automate this process...