#pragma once

#include "TrackMSTargetColor.h"

/**
	RGB color mean-shift target
*/
class trackMSTargetRGB: public TrackMSTargetColor
{
public:	//Member functions

	/**
		Bin assignment function
	*/
	unsigned int quantize(unsigned char r,unsigned  char g,unsigned  char b);

	/**
		Create target candidate from previous candidate
	*/
	void candidate(TrackMSTargetBase*& dest, unsigned char* image, double x, double y, std::vector<double> scale);

	/**
		Updates target candidate without reallocating it
	*/
	void update_candidate(TrackMSTargetBase* dest, unsigned char* image, double x, double y, std::vector<double> scale);

	/**
		Destroys allocated candidate
	*/
	void kill_candidate(TrackMSTargetBase*& dest);

public:	//Constructors/destructor

	trackMSTargetRGB(double x,							//x-coordinate of target ellipse center
					  double y,							//y-coordinate
					  double hx,						//Width of ellipse
					  double hy,						//Height
					  std::vector<double> s,					//Scales to track
					  unsigned char* image,				//Input frame
					  unsigned int imageWidth,			//Width of input frame
					  unsigned int imageHeight,			//Height
					  double backSize = -1,				//Size of background for background reweighting (if -1, no background reweighting is used)
					  unsigned int rNbBin = 16,			//Number of bins in 'red' dimension
					  unsigned int gNbBin = 16,			//Number of bins in 'green' dimension
					  unsigned int bNbBin = 16);		//Number of bins in 'blue' dimension

public: //Operators

};

