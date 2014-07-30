#include <iostream>
#include <vector>

#include "trackFeature.h"
#include "track.h"
#include "trackMatch.h"
#include "trackMeanShift.h"
#include "trackMSTargetRGB.h"

#include <cmath>

//#define MATLAB_INTERFACE		//Matlab interface to capture video frames and display tracked object, comment this if not using matlab
#ifdef MATLAB_INTERFACE
#include "mli.h"
#endif

using namespace std;

int main(int argc, char *argv[], char *envp[])
{
	double x, y;					//x,y coordinates of the target model, those are obtained after manually selecting an object to track
	double hx, hy;					//width and height of the selected object
	unsigned int vsize, hsize;		//Vertical and horizontal sizes of frames
	unsigned char* frame = NULL;	//Buffer for current frame
	unsigned char* roi=NULL;		//Buffer for roi containing the model returned by mli::getTarget
	unsigned char* bw =NULL;		//Non-rectangular model region within roi

	//NOTE: we assume that the frames have INTERLEAVED color channels 

#ifdef MATLAB_INTERFACE
	//Here we use the matlab interface to load movie and get initial frame and target. Replace this with the interface you use
	mli me(atoi(argv[2]),atoi(argv[3]));			//argv[2] is a number specifying which frame to track from, argv[3] specifes how many frames to track for
	me.loadMovie(string(argv[1]));					//argv[1] is the name of the avi file 
	me.getFrame(string(argv[2]),frame,hsize,vsize);	//Get initial frame
	me.getTarget(x,y,hx,hy,roi,bw);					//Here we use the matlab interface to set x,y, hx and hy, roi and bw
#endif

	/*****	0. Setup parameters *****/
	vector<double> s(5);	//Default tracking scales. These should be centered on 1.0. Here we use: 0.8, 0.9, 1.0, 1.1 and 1.2. 
	s[2] = 1;
	s[0] = 0.8;
	s[1] = 0.9;
	s[3] = 1.1;
	s[4] = 1.2;

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
	trackMeanShift::limits bounds(epsSpatial,epsScale,maxAll,maxSpatial,maxScale);		
	
	/***** 1. Create target model histogram *****/
	trackMSTargetRGB model(x,y,hx,hy,s,frame,hsize,vsize,2);

	/***** 1B. Create SIFT feature descriptor for model *****/
	trackSIFTFilter sift(noctaves,nlevels,o_min);		//Declare SIFT Filter object (which contains model features)
	sift.compute(roi,hx,hy,true,bw);					//Create SIFT descriptor for model and prune boundaries

	/***** 1C. Create Matching class for feature descriptors *****/
	trackMatchTri match(mt);

	/***** 2. Create a meanshift tracking object and initialize it with above parameters *****/
	trackMeanShift track1(&model,b,n,bounds,&sift,&match);		
//	trackMeanShift track1(&model,b,n,bounds);

	/***** 3. Create a vector of tracking objects (for now we only use 1 meanshift object but this will be extended) *****/
	vector<trackFlow*> tracker;	
	tracker.push_back(&track1);							//Insert meanshift tracker in vector
 
	/***** 4. Declare output parameters *****/
	pair<double,double> init(0,0);						//This one will contain the x and y coordinates of the target 
	double sc = 0;										//This one will contain the scale "	"	"	"	"	"

#ifdef MATLAB_INTERFACE
	for(int i=0;i<me.nbIter();i++)						//Track for as many iterations as specified by nbIter()
	{
		me.nextFrame(frame);							//Get next frame using matlab interface
#endif

		/***** 5. Track 1 frame: current x,y coordinates and scale are returned in 'init' and 'sc' *****/
		track1Frame(init,sc,frame,tracker,hx,hy,hsize,vsize);			

		tracker[0]->currentCoordinates().show();		//Output current x,y location and scale to stdout (for testing)

#ifdef MATLAB_INTERFACE
		//Display result using matlab interface
		me.showTarget(tracker[0]->currentCoordinates().s_x,tracker[0]->currentCoordinates().s_y,hx,hy,tracker[0]->currentCoordinates().s_scale);
	}
#endif


	//Clean up before exit
	delete[] frame;
	delete[] roi;
	delete[] bw;

	return 0;
}
