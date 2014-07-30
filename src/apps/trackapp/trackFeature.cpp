#include "trackFeature.h"

#include <iostream>

#include <cmath>

extern "C" {
#include <vl/generic.h>
}


void trackFeature::show() const
{
	std::cout<<"Feature point x: "<<m_x<<", y:"<<m_y<<std::endl;
}

void trackFeatArray::show() 
{
	cout<<"Number of feature points: "<<m_feature.size()<<endl;

	for(vector<trackFeature*>::iterator p = m_feature.begin() ; p != m_feature.end(); p++)
	{
		(*p)->show();
	}
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


void trackSIFTFilter::compute(unsigned char* image, unsigned int width, unsigned int height, bool isModel,unsigned char* bw)
{
	m_filter = vl_sift_new(width, height, m_noctaves, m_nlevels, m_nomin);

	unsigned char* gsimage;
	rgb2gray(image,gsimage,width,height);

	//Convert to float array
	vl_sift_pix* pImage = new vl_sift_pix[width*height];
	vl_sift_pix* p = pImage;
	unsigned char* q = gsimage;
	for(unsigned int i=0;i<width*height;i++,p++,q++)
	{
		*p = static_cast<float>(*q);
	}

	vl_sift_process_first_octave(m_filter,pImage);			//Process first octave
	process();

	//Process remaining octaves
	for(int j=0;j<m_noctaves;j++)
	{
		vl_sift_process_next_octave(m_filter);			//Process octave
		process();		
	}

	if (isModel==true)
	{
		//If it is the model, we remove boundary points
		prune(bw,width);
	}

//	process();	//TODO why would there be another call to process if it's already called earlier!


	vl_sift_delete(m_filter);
	m_filter = NULL;

	delete[] pImage;
	delete[] gsimage;
}

void trackSIFTFilter::process()
{
	vl_sift_detect(m_filter);								//Detect keypoints
	int nkp = vl_sift_get_nkeypoints(m_filter) ;				//Get number of keypoints
	m_kp = vl_sift_get_keypoints(m_filter);						//Get the list
	
	double angles[4];

	VlSiftKeypoint * pq;

	unsigned int size = 128;	//TODO: how to determine this empirically?

	for(int i=0;i<nkp;i++)
	{
		int nangle = vl_sift_calc_keypoint_orientations(m_filter,angles,&m_kp[i]);			//Get orientations

		for(int j=0;j<nangle;j++)
		{
			m_pFeature.feature().push_back(new trackFeatureSIFT(m_kp[i].x,m_kp[i].y,size));

			//Here assign descriptor vector of proper length (see docs) which is part of a trackFeatureSIFT object 
			vl_sift_calc_keypoint_descriptor(m_filter,((trackFeatureSIFT*)m_pFeature.feature().back())->desc(),&m_kp[i],angles[j]);
		}
	}
}

void trackSIFTFilter::prune(unsigned char* validity,unsigned int width)
{	

	vector<trackFeature*>& ref = m_pFeature.feature();

	for(vector<trackFeature*>::iterator p = ref.begin();p != ref.end();)
	{
		//Get a feature point
		pair<double,double> xy = ((trackFeatureSIFT*)*p)->spatial();	//TODO: check that this works

		if(validity[(unsigned int)floor(xy.first + xy.second*width)] != 1)
		{
			//This is not a valid point
			//note: this was a reference before, does it still work?
			vector<trackFeature*>::iterator q = p - 1;	//Store current iterator in buffer
			m_pFeature.feature().erase(p);			//Erase current iterator
			p = q;									//Assign previous location
			p++;									//Go to next location
		}
		else
		{
			//This is a valid one
			p++;
		}
	}
}

trackSIFTFilter::trackSIFTFilter(int noctaves, int nlevels, int o_min):m_noctaves(noctaves),m_nlevels(nlevels),m_nomin(o_min)
{
}

trackSIFTFilter::~trackSIFTFilter()
{}
