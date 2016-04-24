#include <iostream>
#include <vector>

#include "TrackFeature.h"
#include "Track.h"
#include "TrackMatch.h"
#include "TrackMeanShift.h"
#include "TrackMSTargetRGB.h"

#include <cmath>

#define OPENCV_INTERFACE 		//To use opencv to load/roi video frames
#ifdef OPENCV_INTERFACE
#include "opencv2/opencv.hpp"
#endif

//#define MATLAB_INTERFACE		//Matlab interface to capture video frames and display tracked object, comment this if not using matlab
#ifdef MATLAB_INTERFACE
#include "mli.h"
#endif

using namespace std;
using namespace cv;

int main(int argc, char *argv[], char *envp[])
{
	double x, y;						//x,y coordinates of the target model, those are obtained after manually selecting an object to track
	double hx, hy;						//width and height of the selected object
	int vsize, hsize;					//Vertical and horizontal sizes of frames
	unsigned char* roi = nullptr;		//Buffer for roi containing the model returned by mli::getTarget
	//unsigned char* bw = nullptr;		//Non-rectangular model region within roi
	vector<unsigned char> bw;

#ifdef MATLAB_INTERFACE
	unsigned char* frame = nullptr;	//Buffer for current frame
	//Here we use the matlab interface to load movie and get initial frame and target. Replace this with the interface you use
	mli me(atoi(argv[2]),atoi(argv[3]));			//argv[2] is a number specifying which frame to track from, argv[3] specifes how many frames to track for
	me.loadMovie(string(argv[1]));					//argv[1] is the name of the avi file 
	me.getFrame(string(argv[2]),frame,hsize,vsize);	//Get initial frame
	me.getTarget(x,y,hx,hy,roi,bw);					//Here we use the matlab interface to set x,y, hx and hy, roi and bw
#endif

#ifdef OPENCV_INTERFACE
	string movie_path = argv[1];
	VideoCapture vid(movie_path);

	namedWindow("image", 1);
	namedWindow("target", 2);
	Mat frame;
	vid >> frame;
	Size size(640, 480);
	Mat resized_frame;
	resize(frame, resized_frame, size);
	imshow("image", resized_frame);

	//Extract rectangle
	x = y = 250;
	hx = hy = 100;
	vsize = resized_frame.rows;
	hsize = resized_frame.cols;
	Mat target = resized_frame(Range(x - hx/2, x + hx/2), Range(y - hy/2, y + hy/2));
	imshow("target", target);
	waitKey();

	bw.resize(hx * hy);
	memcpy(&bw[0], target.data, hx * hy);
#endif

	vector<double> s = {0.8, 0.9, 1.0, 1.1, 1.2};	//Default tracking scales. These should be centered on 1.0. Here we use: 0.8, 0.9, 1.0, 1.1 and 1.2.

	// More parameters to be set by that user (here use default values)
	double b            = 1.1;				
	unsigned int n      = 2;
	unsigned int maxAll = 20, maxSpatial=15, maxScale=15;
	double epsSpatial   = 0.2;
	double epsScale     = 0.1;
	
	// Parameters for SIFT descriptors
	int noctaves		= 2;	//Number of octaves
	int nlevels		    = 3;	//Number of levels per octave
	int o_min			= 1;	//Minimum number of orientations

	//Parameters for matching
	double mt			= 0.1;	//Threshold at which the first neighbor must be higher than the second one

	//Create a 'limits' object: this governs convergence properties and can be used as parameter to determine how good the quality of the solution is
	TrackMeanShift::limits bounds(epsSpatial, epsScale, maxAll, maxSpatial, maxScale);

	//Create target histogram
	TrackMSTargetRGB model(x, y, hx, hy, s, frame.data, hsize, vsize, 2);

	//Create SIFT descriptor
	TrackSIFTFilter sift(noctaves, nlevels, o_min);	//Declare SIFT Filter object (which contains model features)
	sift.compute(target.data, hx, hy, true, bw); //Create SIFT descriptor for model and prune boundaries

	//Create Matching class for feature descriptors
	trackMatchTri match(mt);

	//Create a meanshift tracking object and initialize it with above parameters
	int unique_id = 0;
	TrackMeanShift track1(unique_id, &sift, &match, &model, b, n, bounds);

	//Create a vector of tracking objects (for now we only use 1 meanshift object but this will be extended)
	vector<TrackFlow*> tracker;
	tracker.push_back(&track1);							//Insert meanshift tracker in vector

	//Declare output parameters
	auto init = make_pair(0.0, 0.0);	//This one will contain the x and y coordinates of the target
	double sc = 0;	//This one will contain the scale "	"	"	"	"	"

#ifdef MATLAB_INTERFACE
	for(int i = 0;i < me.nbIter(); i++)						//Track for as many iterations as specified by nbIter()
	{
		me.nextFrame(frame);							//Get next frame using matlab interface
#endif

#ifdef OPENCV_INTERFACE
	namedWindow("image", 1);
	for(int i = 0; i < /*vid.get(cv::CAP_PROP_FRAME_COUNT)*/ 10; i++)						//Track for as many iterations as specified by nbIter()
	{
		cout<<"Frame "<<i<<endl;
		//Mat frame;
		vid >> frame;
		//Size size(640, 480);
		//Mat resized_frame;
		resize(frame, resized_frame, size);
		imshow("image", resized_frame);
		waitKey(30);
#endif

		// Track 1 frame: current x,y coordinates and scale are returned in 'init' and 'sc'
		track_one_frame(init, sc, frame.data, tracker, hx, hy, hsize, vsize);
//
//		tracker[0]->currentCoordinates().show();		//Output current x,y location and scale to stdout (for testing)


#ifdef MATLAB_INTERFACE
		//Display result using matlab interface
		me.showTarget(tracker[0]->currentCoordinates().s_x,tracker[0]->currentCoordinates().s_y,hx,hy,tracker[0]->currentCoordinates().s_scale);
	}
#endif

#ifdef OPENCV_INTERFACE
	}
#endif

	//Clean up before exit

#ifdef MATLAB_INTERFACE
	delete[] frame;
	delete[] bw;
#endif
	delete[] roi;

#ifdef OPENCV_INTERFACE
	vid.release();
#endif

	return 0;
}
