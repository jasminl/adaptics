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

void TrackFeatArray::show() 
{
	cout<<"Number of feature points: "<<_feature.size()<<endl;

	for(vector<TrackFeature*>::iterator p = _feature.begin() ; p != _feature.end(); p++)
	{
		(*p)->show();
	}
}


double TrackFeatureSIFT::compare(TrackFeature* x)
{
	TrackFeatureSIFT* px = dynamic_cast<TrackFeatureSIFT*>(x);
	double d = 0; //todo check if gets caught by vectorizer
	for(int i = 0; i < (int)_desc.size(); i++)
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
		*pout = static_cast<unsigned char>(.3 * *pimage + .59 * *(pimage+1) + .11 * *(pimage+2));	
}


void TrackSIFTFilter::compute(unsigned char* image, unsigned int width, unsigned int height,
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

void TrackSIFTFilter::process()
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
			_p_feat.feature().push_back(new TrackFeatureSIFT(_kp[i].x,_kp[i].y,size));

			//Here assign descriptor vector of proper length (see docs) which is part of a trackFeatureSIFT object 
			vl_sift_calc_keypoint_descriptor(_filter,((TrackFeatureSIFT*)_p_feat.feature().back())->desc(),&_kp[i],angles[j]);
		}
	}
}

void TrackSIFTFilter::prune(unsigned char* validity, unsigned int width)
{	
	auto& ref = _p_feat.feature();
	cout<<_p_feat.size()<<endl;
	for(auto p = ref.begin();p != ref.end();)
	{
		//Get a feature point

		auto xy = ((TrackFeatureSIFT*)*p)->spatial();	//TODO: check that this works
		cout<<xy.first<<" "<<xy.second;

		if(validity[(int)floor(xy.first + xy.second * width)] != 1)
		{	//This is not a valid point
			//note: this was a reference before, does it still work?
			cout<<" invalid "<<endl;
			vector<TrackFeature*>::iterator q = p - 1;	//Store current iterator in buffer
			_p_feat.feature().erase(p);			//Erase current iterator
			p = q;									//Assign previous location
			p++;									//Go to next location
		}
		else
		{   //This is a valid one
			cout<<" valid"<<endl;
			p++;
		}
	}
}

TrackSIFTFilter::TrackSIFTFilter(int noctaves, int nlevels, int o_min)
: _filter(nullptr), _kp(nullptr), _noctaves(noctaves), _nlevels(nlevels),
  _nomin(o_min)
{}

TrackSIFTFilter::~TrackSIFTFilter()
{}
