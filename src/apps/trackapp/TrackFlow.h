#pragma once

#include <vector>

#include "TrackFeature.h"
#include "TrackMatch.h"

/**
	Base class for all tracking methods
*/
class TrackFlow
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
			Constructor
		*/
		_coord(double x, double y, double scale, std::vector<double> scaleRange)
		: s_x(x), s_y(y), s_scale(scale), s_scaleRange(scaleRange)
		{}

		_coord(): s_x(0.0), s_y(0.0), s_scale(0.0), s_scaleRange(std::vector<double>())
		{}

		/**
			Display parameters
		*/
		void show() const;

	}coord;

public:

	/**
	 * Base constructor
	 */
	TrackFlow(TrackFilter* filt, TrackMatch* match)
	: _cur_x(0), _cur_y(0), _cur_scale(0), _feat_filter(filt),
	  _target_filter(filt), _match(match)
	{}

public:

	virtual
	bool track(unsigned char* image);

	coord currentCoordinates() const;

	TrackFilter* featFilter();

	TrackFilter*& targetFilter();

	TrackMatch* matcher();

protected:

	/*<< The following variables are returned by tracking methods */
	double _cur_x;							/*<< Current x-coordinates of center of target (for some methods only) */
	double _cur_y;							/*<<"	"	 y-"	"	"	"	"	"	"	"	"	"	"	"	"	"  */
	double _cur_scale;						/*<<:	:	scale	"	"	"	"	"	"	"	"	"	"	"	"	"  */

	TrackFilter*   _feat_filter;				/*<< Feature filter for the model */
	TrackFilter*   _target_filter;				/*<< Feature filter for a particular target */
	TrackMatch*    _match;						/*<<Matching object */

	std::vector<double> _cur_expanded_scale;		/*<< Current range of scales tracked */
	coord _coord;								/*<< To return coordinates */
};

inline
bool TrackFlow::track(unsigned char* image)
{
	return true;
}

inline
TrackFlow::coord TrackFlow::currentCoordinates() const
{
	return _coord;
}

inline
TrackFilter* TrackFlow::featFilter()
{
	return _feat_filter;
}

inline
TrackFilter*& TrackFlow::targetFilter()
{
	return _target_filter;
}

inline
TrackMatch* TrackFlow::matcher()
{
	return _match;
}
