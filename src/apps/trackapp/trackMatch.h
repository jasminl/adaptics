#pragma once

#include <vector>

#ifdef TRACK_EXPORTS
#define TRACK_API __declspec(dllexport)
#else
#define TRACK_API __declspec(dllimport)
#endif

using namespace std;

class trackFeatArray;
class trackFeature;

/**
	Base class for matching point pairs
*/
class  TRACK_API trackMatch
{
public:

	/**
		Default constructor
	*/	
	trackMatch(void){}

	/**
		Destructor
	*/
	virtual ~trackMatch(void){}

	/**
		Main matching function
	*/
	virtual void match(trackFeatArray* model, trackFeatArray* target)=0;

protected:	//Data members

	//Matched points are aligned in those two vectors which store x,y coordinates
	vector< pair<double,double> > m_x1;	/*<<Model points */
	vector< pair<double,double> > m_x2; /*<<Target points */

	double m_rotation;	/*<<Rotation angle */
	double m_xt;		/*<<Horizontal translation */
	double m_yt;		/*<<Vertical translation */
	double m_xh;		/*<<Horizontal scale change */
	double m_yh;		/*<<Vertical scale change */

};

/**
	Simple matching
	See: Chen et al. (2008) Mean-shift tracking combining SIFT. ICSP2008.
*/
class  TRACK_API trackMatchTri: public trackMatch
{
public:	//Member functions

	/**
		Main matching function
	*/
	void match(trackFeatArray* model, trackFeatArray* target);

	/**
		Find 1st and 2nd nearest neighbors
	*/
	pair< pair<unsigned int, double>, pair<unsigned int, double> > firstSecond(trackFeature* modelPt, trackFeatArray* target);

public:	//Constructors / destructor

	/**
		Default constructor
	*/
	trackMatchTri(){}

	/**
		Constructor with threshold parameter
	*/
	trackMatchTri(double threshold):trackMatch(),m_t(threshold)
	{}

	/**
		Destructor
	*/
	virtual ~trackMatchTri(){}

protected:	//Data members

	double m_t;		/*<<Threshold for ratio of first to second neighbor difference */
};