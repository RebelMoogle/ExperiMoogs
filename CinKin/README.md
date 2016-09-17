# Cinder with Kinect V2 SDK
-------------4

This is an experiment to get Kinect working with Cinder V2

## Requisites: 
* Windows 8/10 ONLY :(
* Visual Studio 2013 
    * VS2015 works too, but you need to get msvcr120d.dll and msvcd120.dll yourself (can't legally share, it, Google is your friend) 
* [Kinect SDK V2](https://www.microsoft.com/en-us/download/details.aspx?id=44561)
* [Kinect V2 Common Bridge Cinder Block](https://github.com/wieden-kennedy/Cinder-KCB2)
    * [VS2015 Version](https://github.com/wieden-kennedy/Cinder-KCB2/tree/vc2015)
* [Cinder](https://libcinder.org/download)

1. Install Cinder, VS2015 and Kinect SDK
2. Copy the Kinect Common Bridge block to CINDER_DIR\blocks\ (in it's own folder)
3. Start Tinderbox and create a VS 2013 project 
    * Select KCB2: Basic Template
4. (VS2015 only) Copy msvcr120d.dll to your windows\system32 folder
5. Build and run! (make sure your kinect is connected and installed.)