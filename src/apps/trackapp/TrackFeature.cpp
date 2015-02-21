#include "TrackFeature.h"

#include <iostream>

#include <cmath>

extern "C" {
#include <vl/generic.h>
}

using namespace std;

void TrackFeature::show() const
{
	cout<<"Feature point x: "<<_x<<", y:"<<_y<<endl;
}

void trackFeatArray::show() 
{
	cout<<"Number of feature points: "<<_feature.size()<<endl;

	for(vector<TrackFeature*>::iterator p = _feature.begin() ; p != _feature.end(); p++)
	{
		(*p)->show();
	}
}


double trackFeatureSIFT::compare(TrackFeature* x)
{
	trackFeatureSIFT* px = dynamic_cast<trackFeatureSIFT*>(x);
	double d = 0;
	for(unsigned int i = 0; i < _size; i++)
		d += (_desc[i] - px->desc()[i]) * (_desc[i] - px->desc()[i]);
	return d;
}

/**
	Converts from rgb to grayscale, assumes interleaved channels
*/
void rgb2gray(unsigned char* image, unsigned char*& output, unsigned int width, unsigned int height)
{
	output = new unsigned char[width*height];

	unsigned char* pimage = image;
	unsigned char* pout = output;

	for(unsigned int i=0;i<width*height;i++, pout++, pimage+=3)
	{
		*pout = static_cast<unsigned char>(.3 * *pimage + .59 * *(pimage+1) + .11 * *(pimage+2));	
	}
}


void trackSIFTFilter::compute(unsigned char* image, unsigned int width, unsigned int height,
		bool isModel,unsigned char* bw)
{
	_filter = vl_sift_new(width, height, _noctaves, _nlevels, _nomin);

	unsigned char* gsimage;
	rgb2gray(image, gsimage, width, height);

	//Convert to float array
	vl_sift_pix* pImage = new vl_sift_pix[width * height];
	vl_sift_pix* p = pImage;
	unsigned char* q = gsimage;
	for(unsigned int i = 0; i < width * height; i++, p++, q++)
		*p = static_cast<float>(*q);

	vl_sift_process_first_octave(_filter, pImage); //Process first octave
	process();

	//Process remaining octaves
	for(int j = 0; j < _noctaves; j++)
	{
		vl_sift_process_next_octave(_filter); //Process octave
		process();		
	}

	if (isModel == true)//If it is the model, we remove boundary points
		prune(bw, width);

//	process();	//TODO why would there be another call to process if it's already called earlier!


	vl_sift_delete(_filter);
	_filter = nullptr;

	delete[] pImage;
	delete[] gsimage;
}

void trackSIFTFilter::process()
{
	vl_sift_detect(_filter);								//Detect keypoints
	int nkp = vl_sift_get_nkeypoints(_filter) ;				//Get number of keypoints
	_kp = vl_sift_get_keypoints(_filter);						//Get the list
	
	double angles[4];
	unsigned int size = 128;	//TODO: how to determine this empirically?

	for(int i=0;i<nkp;i++)
	{
		int nangle = vl_sift_calc_keypoint_orientations(_filter,angles,&_kp[i]);			//Get orientations

		for(int j=0;j<nangle;j++)
		{
			m_pFeature.feature().push_back(new trackFeatureSIFT(_kp[i].x,_kp[i].y,size));

			//Here assign descriptor vector of proper length (see docs) which is part of a trackFeatureSIFT object 
			vl_sift_calc_keypoint_descriptor(_filter,((trackFeatureSIFT*)m_pFeature.feature().back())->desc(),&_kp[i],angles[j]);
		}
	}
}

void trackSIFTFilter::prune(unsigned char* validity,unsigned int width)
{	

	vector<TrackFeature*>& ref = m_pFeature.feature();

	for(vector<TrackFeature*>::iterator p = ref.begin();p != ref.end();)
	{
		//Get a feature point
		auto xy = ((trackFeatureSIFT*)*p)->spatial();	//TODO: check that this works

		if(validity[(unsigned int)floor(xy.first + xy.second*width)] != 1)
		{
			//This is not a valid point
			//note: this was a reference before, does it still work?
			vector<TrackFeature*>::iterator q = p - 1;	//Store current iterator in buffer
			m_pFeature.feature().erase(p);			//Erase current iterator
			p = q;									//Assign previous location
			p++;									//Go to next location
		}
		else
		{   //This is a valid one
			p++;
		}
	}
}

trackSIFTFilter::trackSIFTFilter(int noctaves, int nlevels, int o_min)
: _filter(nullptr), _kp(nullptr), _noctaves(noctaves), _nlevels(nlevels),
  _nomin(o_min)
{}

trackSIFTFilter::~trackSIFTFilter()
{}
