#pragma once

#include <vector>
#include "TrackFlow.h"

/**
	Function to track 1 frame
*/
 std::vector<bool> track_one_frame(std::pair<double,double>& location,			//Output location
						   double& scale,										//Output scale
						   unsigned char* image,								//image frame
						   std::vector<TrackFlow*> tracker,							//trackFlow objects
						   double hx,											//Width of model
						   double hy,											//Height of model
						   unsigned int width,									//Width of image
						   unsigned int height);								//Height of image


/**
	Function to crop image to region of interest, for feature matching (ASSUMES INTERLEAVED RGB)
*/
unsigned char* crop(unsigned char*& buffer,							//Buffer to store cropped result
							  unsigned int width,								//Desired width of crop
							  unsigned int height,								//Desired height "	"
							  double x,											//x-coordinate of center of part to crop
							  double y,											//y-coordinate of center of "	"	"
							  unsigned char* image,								//Input image from where to crop
							  unsigned int imageWidth,							//Image width
							  unsigned int imageHeight);						//Image height
