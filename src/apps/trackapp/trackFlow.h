#pragma once

#include <vector>

#include "trackFeature.h"
#include "TrackMatch.h"

/**
	Base class for all tracking methods
*/
class  trackFlow
{
public: //Typedefs

	/**
		Structure used to retrieve parameters
	*/
	typedef struct _coord
	{
		double s_x;
		double s_y;
		double s_scale;
		std::vector<double> s_scaleRange;

		/**
			Default constructor
		*/
		_coord(){}

		/**
			Constructor
		*/
		_coord(double x, double y, double scale, std::vector<double> scaleRange):s_x(x),s_y(y),s_scale(scale),s_scaleRange(scaleRange)
		{}

		/**
			Display parameters
		*/
		void show() const;

	}coord;

public:	//Member functions

	virtual bool track(unsigned char* image)
	{
		return true;
	}

	coord currentCoordinates() const
	{
		return m_coord;
	}

	trackFilter* featFilter()
	{
		return m_featFilter;
	}

	trackFilter*& targetFilter()
	{
		return m_targetFilter;
	}

	TrackMatch* matcher()
	{
		return m_match;
	}

public:	//Constructors/destructors
	trackFlow(trackFilter* filt, TrackMatch* match);

protected:	//Data members

	/*<< The following variables are returned by tracking methods */
	double m_currentX;							/*<< Current x-coordinates of center of target (for some methods only) */
	double m_currentY;							/*<<"	"	 y-"	"	"	"	"	"	"	"	"	"	"	"	"	"  */
	double m_currentScale;						/*<<:	:	scale	"	"	"	"	"	"	"	"	"	"	"	"	"  */

	trackFilter*   m_featFilter;				/*<< Feature filter for the model */
	trackFilter*   m_targetFilter;				/*<< Feature filter for a particular target */
	TrackMatch*    m_match;						/*<<Matching object */

	std::vector<double> m_currentExpandedScale;		/*<< Current range of scales tracked */

	coord m_coord;								/*<< To return coordinates */

};
