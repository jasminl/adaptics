#pragma once

#include "TrackMatch.h"

/**
	Base class for random consensus matching algorithms
*/
class   TrackSAC
{
public:

	/**
		Defines rotation, scaling and translation transforms
	*/
	typedef struct _RST
	{
		double _rot;	//!< Rotation angle
		double _scale;	//!< Scaling
		double _t1;		//!< Horizontal translation
		double _t2;		//!< Vertical translation

		/**
			Default constructor, sets all member variables to zero
		*/
		_RST()
		: _rot(0.0), _scale(0.0), _t1(0.0), _t2(0.0)
		{}

		/**
			Constructor
			\param rot Rotation angle
			\param scale Transform scale
			\param t1 Horizontal translation
			\param t2 Vertical translation
		*/
		_RST(double rot, double scale, double t1, double t2)
		:_rot(rot), _scale(scale), _t1(t1), _t2(t2)
		{}

	}rst;

public:

	/**
		Default constructor
	*/
	TrackSAC(void)
	{
		_match = nullptr;
		_alarm_rate = 0.000001;
		_min_iter = 100;			
		_max_iter = 100000;
		_adapt_iter = 0;
	}

	/**
		Constructor with matched point object
		\param pt Matching point
	*/
	TrackSAC(TrackMatch* pt, double alarm_rate = 0.000001, unsigned int min_iter = 100,
		unsigned int max_iter = 10000)
	:_match(pt), _alarm_rate(alarm_rate),_min_iter(min_iter),_max_iter(max_iter),
	 _adapt_iter(0)
	{}

	/**
		Destructor, performs no deallocation.
	*/
	virtual
	~TrackSAC(void)
	{}

public:

	/**
		Estimates parameters for rotation, scale and translation
	*/
	virtual
	void estimate_RST(rst& out){}	//TODO: maybe this shouldn't be derived...
	
	/**
		Estimates number of iterations to be run
	*/
	double estimate_iter(double probQ);

	/**
		Estimates inlier MSS probability
	*/
	double estimate_inlier_MSS(unsigned int N, unsigned int NInlier, unsigned int k);

	/**
		Get minimal sample set
	*/
	void mss(rst& out){};//TODO: complete this

protected:
	
	TrackMatch*	_match;			//!< Matched points
	double _alarm_rate;				//!< Threshold to determine number of iterations
	unsigned int _min_iter;			//!< Minimum number of iterations to perform SAC
	unsigned int _max_iter;			//!< Maximum number of iterations to perform SAC
	unsigned int _adapt_iter;		//!< Number of iterations to perform based on calculations
};

class   TrackRANSAC: public TrackSAC
{

public:
	
	/**
		Default constructor
	*/
	TrackRANSAC(){}

	/**
		Constructor with matched points
	*/
	TrackRANSAC(TrackMatch* pt, double alarmRate)
	: TrackSAC(pt, alarmRate)
	{}

public:
	void estimate_RST(rst& out);
};
