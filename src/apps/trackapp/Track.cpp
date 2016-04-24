#include <numeric>
#include <algorithm>
#include <iostream>
#include <string.h>

#include "Track.h"
#include "TrackFeature.h"
#include "TrackMatch.h"

using namespace std;


vector<bool> track_one_frame(pair<double,double>& location, double& scale, unsigned char* image,
		vector<TrackFlow*> tracker, double hx, double hy, unsigned int width, unsigned int height)
{
	//Create boolean vector indicating whether the tracked object is valid, for each tracking object
 	vector<bool> validTrack(tracker.size());
	auto q = validTrack.begin();

	double denominator = 0;				//Qties used to average output
	location.first = location.second = scale = 0;	//Initialize output variables to zero

	for(auto p: tracker)
	{
		*q = p->track(image);

		if(*q == true)
		{
			//If this a valid target, include it in the fused x, y, scale estimates.
			TrackFlow::coord current = p->currentCoordinates();
			location.first += current.s_x;
			location.second += current.s_y;
			scale += current.s_scale;

			denominator++;
		}
	}

	//Average output
	location.first /= denominator;
	location.second /= denominator;
	scale /= denominator;

	//Create cropped image for shape matching
	unsigned char* cropped = nullptr;
	
	crop(cropped, scale*hx,scale*hy,location.first,location.second,image,width,height);
	
	//Perform shape matching if applicable
	for(auto p = tracker.begin() ; p != tracker.end() ; p++)
	{
		if((*p)->featFilter() == nullptr) continue;		//Skip this step if no feature descriptor has been defined

		//Use target filter to compute and process keypoints
		vector<unsigned char> empty; //todo remove this..
		(*p)->targetFilter()->compute(cropped, scale * hx, scale * hy, false, empty);
	//	(*p)->targetFilter()->process();
		
		//DEBUGGIN
		(*p)->featFilter()->feature()->show();		//THIS SHOWS that zero points are in the array
		(*p)->targetFilter()->feature()->show();

//TODO: show features on image
		//Perform matching between target and model
		(*p)->matcher()->match((*p)->featFilter()->feature(),(*p)->targetFilter()->feature());

//TODO: show matches on image
	}

	delete[] cropped;	//Clean up matching image

	return validTrack;	//Return vector indicating validity of tracked object
}

unsigned char* crop(unsigned char*& buffer, unsigned int width, unsigned int height,double x,double y,unsigned char* image,unsigned int imageWidth,unsigned int imageHeight)
{
	if (buffer != nullptr)
	{
		delete[] buffer;
	}

	//Determine coordinates of area to crop
	unsigned int x1 = max((double)0,x-width/2);
	unsigned int x2 = min((double)imageWidth,x+width/2);
	unsigned int y1 = max((double)0,y-height/2);
	unsigned int y2 = min((double)imageHeight,y+height/2);

	//Calculate range
	double xRange = x2 - x1 + 1;
	double yRange = y2 - y1 + 1;

	//Allocate buffer
	buffer = new unsigned char[3 * (unsigned int)(ceil(xRange) * ceil(yRange))];

	//Copy image
	unsigned char* p = buffer;								//Point to cropped image
	unsigned char* q = image + y1 * imageWidth + x1;		//Pointer to base image
	for(unsigned int i=y1;i<y2;i++)
	{
		memcpy(p,q,3 * xRange * sizeof(unsigned char));		//Copy 1 row
		q += 3 * imageWidth;								//Forward to next row
	}

	return buffer;											//Return result
}


void TrackFlow::_coord::show() const
{
	cout<<"(X,Y)=("<<s_x<<","<<s_y<<") Scale="<<s_scale<<" Range=(";

	for(unsigned int i=0;i<s_scaleRange.size(); i++) cout<<s_scaleRange[i]<<" ";

	cout<<")\n";
}
