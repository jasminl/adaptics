#include "TrackMSTargetRGB.h"

#include <cmath>

using namespace std;

unsigned int TrackMSTargetRGB::quantize(unsigned char r,unsigned  char g,unsigned  char b)
{
	unsigned int red   = static_cast<unsigned int>(floor(static_cast<double>(r)/_levels_per_bin[0]));
	unsigned int green = static_cast<unsigned int>(floor(static_cast<double>(g)/_levels_per_bin[1]));
	unsigned int blue  = static_cast<unsigned int>(floor(static_cast<double>(b)/_levels_per_bin[2]));

	return red + green * _bins_per_dim[0] + blue * _bins_per_dim[0]*_bins_per_dim[1];
}

void TrackMSTargetRGB::candidate(TrackMSTargetBase*& dest, unsigned char* image, double x, double y, std::vector<double> scale)
{
	kill_candidate(dest);
	dest = new TrackMSTargetRGB(x,y,_hx,_hy,scale,image,_width,_height,_back_size,_bins_per_dim[0],_bins_per_dim[1],_bins_per_dim[2]);
}


/***** trackMSTargetRGB *****/

TrackMSTargetRGB::TrackMSTargetRGB(double x,
					  double y,
					  double hx,
					  double hy,
					  vector<double> s,
					  unsigned char* image,
					  unsigned int imageWidth,
					  unsigned int imageHeight,
					  double backSize,
					  unsigned int rNbBin,
					  unsigned int gNbBin,
					  unsigned int bNbBin):TrackMSTargetColor(x,y,hx,hy,s,image,imageWidth,imageHeight,backSize)
{

	_bins_per_dim.push_back(rNbBin);										//Setup histogram parameters
	_bins_per_dim.push_back(gNbBin);
	_bins_per_dim.push_back(bNbBin);

	_levels_per_bin.push_back(256.0 / _bins_per_dim[0]);				//Get number of color levels per bin for each dimension
	_levels_per_bin.push_back(256.0 / _bins_per_dim[1]);
	_levels_per_bin.push_back(256.0 / _bins_per_dim[2]);

	allocateHistograms(s,(backSize>0)?true:false);					//Allocate histogram space

	makeHistogram(_s);												//Make histograms
}


void TrackMSTargetRGB::update_candidate(TrackMSTargetBase* dest, unsigned char* image, double x, double y, std::vector<double> scale)
{
	TrackMSTargetRGB* p = dynamic_cast<TrackMSTargetRGB*>(dest);

	std::pair< std::vector<unsigned int> , std::vector<double> > q = p->histParams();

	if(_bins_per_dim != q.first)				//Check that histogram doesn't change
	{
		p->set_histogram(_bins_per_dim,_levels_per_bin);
	}

	p->allocateHistograms(scale,(_back_size>0)?true:false);	//Always reallocate to keep track of scale range change

	p->allocate(scale);							//Allocate remaining base histograms

	p->set(x,y,scale);							//Set new space scale location

	p->makeHistogram(scale);					//Make histogram
}

void TrackMSTargetRGB::kill_candidate(TrackMSTargetBase*& dest)
{
	delete dynamic_cast<TrackMSTargetRGB*>(dest);
	dest = NULL;
}
