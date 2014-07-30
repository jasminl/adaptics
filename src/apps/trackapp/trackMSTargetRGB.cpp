#include "trackMSTargetRGB.h"

#include <cmath>

using namespace std;

unsigned int trackMSTargetRGB::quantize(unsigned char r,unsigned  char g,unsigned  char b)
{
	unsigned int red   = static_cast<unsigned int>(floor(static_cast<double>(r)/m_nbLevelsPerBin[0]));
	unsigned int green = static_cast<unsigned int>(floor(static_cast<double>(g)/m_nbLevelsPerBin[1]));
	unsigned int blue  = static_cast<unsigned int>(floor(static_cast<double>(b)/m_nbLevelsPerBin[2]));

	return red + green * m_nbBinsPerDim[0] + blue * m_nbBinsPerDim[0]*m_nbBinsPerDim[1];
}

void trackMSTargetRGB::candidate(TrackMSTargetBase*& dest, unsigned char* image, double x, double y, std::vector<double> scale)
{
	kill_candidate(dest);
	dest = new trackMSTargetRGB(x,y,_hx,_hy,scale,image,_width,_height,m_backSize,m_nbBinsPerDim[0],m_nbBinsPerDim[1],m_nbBinsPerDim[2]);
}


/***** trackMSTargetRGB *****/

trackMSTargetRGB::trackMSTargetRGB(double x,
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
					  unsigned int bNbBin):trackMSTargetColor(x,y,hx,hy,s,image,imageWidth,imageHeight,backSize)
{

	m_nbBinsPerDim.push_back(rNbBin);										//Setup histogram parameters
	m_nbBinsPerDim.push_back(gNbBin);
	m_nbBinsPerDim.push_back(bNbBin);

	m_nbLevelsPerBin.push_back(256.0 / m_nbBinsPerDim[0]);				//Get number of color levels per bin for each dimension
	m_nbLevelsPerBin.push_back(256.0 / m_nbBinsPerDim[1]);
	m_nbLevelsPerBin.push_back(256.0 / m_nbBinsPerDim[2]);

	allocateHistograms(s,(backSize>0)?true:false);					//Allocate histogram space

	makeHistogram(_s);												//Make histograms
}


void trackMSTargetRGB::update_candidate(TrackMSTargetBase* dest, unsigned char* image, double x, double y, std::vector<double> scale)
{
	trackMSTargetRGB* p = dynamic_cast<trackMSTargetRGB*>(dest);

	std::pair< std::vector<unsigned int> , std::vector<double> > q = p->histParams();

	if(m_nbBinsPerDim != q.first)				//Check that histogram doesn't change
	{
		p->setHistogram(m_nbBinsPerDim,m_nbLevelsPerBin);
	}

	p->allocateHistograms(scale,(m_backSize>0)?true:false);	//Always reallocate to keep track of scale range change

	p->allocate(scale);							//Allocate remaining base histograms

	p->set(x,y,scale);							//Set new space scale location

	p->makeHistogram(scale);					//Make histogram
}

void trackMSTargetRGB::kill_candidate(TrackMSTargetBase*& dest)
{
	delete dynamic_cast<trackMSTargetRGB*>(dest);
	dest = NULL;
}
