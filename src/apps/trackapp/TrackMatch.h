#pragma once

#include <vector>
#include "TrackFeature.h"

/**
	Base class for matching point pairs
*/
class TrackMatch
{
public:

	/**
		Default constructor, perform no allocation
	*/	
	TrackMatch()
	{}

	/**
		Destructor, deallocates no memory
	*/
	virtual
	~TrackMatch()
	{}

	/**
		Main matching function
		\param model Candidate matching model
		\param target Target matched to the model
	*/
	virtual
	void match(trackFeatArray* model, trackFeatArray* target) = 0;

protected:

	//Matched points are aligned in those two vectors which store x,y coordinates
	std::vector<std::pair<double, double>> _x1;	//!< Model points
	std::vector<std::pair<double, double>> _x2; //!< Target points

	double _rotation;	//!< Rotation angle
	double _xt;		//!< Horizontal translation
	double _yt;		//!< Vertical translation
	double _xh;		//!< Horizontal scale change
	double _yh;		//!< Vertical scale change
};

/**
	Simple matching
	See: Chen et al. (2008) Mean-shift tracking combining SIFT. ICSP2008.
*/
class   trackMatchTri: public TrackMatch
{
public:

	/**
		Main matching function
	*/
	void match(trackFeatArray* model, trackFeatArray* target);

	/**
		Find 1st and 2nd nearest neighbors
	*/
	std::pair<std::pair<unsigned int, double>, std::pair<unsigned int, double> > firstSecond(TrackFeature* modelPt, trackFeatArray* target);

public:

	/**
		Default constructor
	*/
	trackMatchTri(){}

	/**
		Constructor with threshold parameter
	*/
	trackMatchTri(double threshold)
	: _t(threshold)
	{}

	/**
		Destructor
	*/
	virtual
	~trackMatchTri()
	{}

protected:
	double _t;		//!< Threshold for ratio of first to second neighbor difference
};
