#pragma once

#include <vector>
#include <map>
#include "fastmap.h"
#include "TrackMSTargetBase.h"

/**
	Class for representing a mean-shift COLOR histogram target
*/
class trackMSTargetColor: public TrackMSTargetBase
{
protected:	//Data members

	fastmap<double,double> m_hist;

	std::map<double, std::vector<double> > m_histogram;	/**< Color histogram, one per scale */
	std::vector<double> m_backHistogram;				/**< Background color histogram (optional) */

	std::vector<unsigned int> m_nbBinsPerDim;		/**< The number of bins per color dimension */
	std::vector<double> m_nbLevelsPerBin;			/**< The corresponding number of color levels in a single bin, per dimension */

	double m_backSize;							/**< Areal factor that determines the size of the background area (relative to target) to use for background weighting */

public:	//Member functions

	/**
		Returns histogram parameters
	*/
	std::pair< std::vector<unsigned int> , std::vector<double> > histParams() const;

	/**
		Weight for target pixels
	*/
	void weight(TrackMSTargetBase* model);

	/**
		Metric based on Bhattacharyya coefficient
	*/
	double bhattacharyya_distance(TrackMSTargetBase* a, double scale);

	/**
		Bin assignment function. For color models that don't have 3 components, set the third one to zero
	*/
	virtual unsigned int quantize(unsigned char q1,unsigned char q2,unsigned char q3)=0;

	/**
		Allocate histogram space
	*/
	void allocateHistograms(std::vector<double> scale, bool backgroundWeight);

	/**
		Histogram building function, assumes interleaved image
	*/
	void makeHistogram(std::vector<double>& scale);


	/**
		Create target candidate from previous candidate
	*/
	virtual void candidate(TrackMSTargetBase*& dest, unsigned char* image, double x, double y, std::vector<double> scale)=0;

	/**
		Updates target candidate
	*/
	virtual void update_candidate(TrackMSTargetBase* dest, unsigned char* image, double x, double y, std::vector<double> scale)=0;

	/**
		Destroy candidate
	*/
	virtual void kill_candidate(TrackMSTargetBase*& dest)=0;

	/**
		Sets histogram parameters
	*/
	void setHistogram(std::vector<unsigned int> nb, std::vector<double> sz);

public: //Constructors/destructor

	/**
		Default constructor for color target representation
	*/
	trackMSTargetColor(double x,
					  double y,
					  double hx,
					  double hy,
					  std::vector<double> s,
					  unsigned char* image,
					  unsigned int imageWidth,
					  unsigned int imageHeight,
					  double backSize);
};


