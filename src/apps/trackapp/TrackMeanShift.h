#pragma once

#include "TrackMSTargetBase.h"
#include "TrackFlow.h"

/**
	Main class for meanshift tracking
*/
class  TrackMeanShift: public TrackFlow
{
public: //typedefs etc.

	/**
		All limit constants for meanshift tracking in scale-space (note, scale tracking is optional)
	*/
	typedef struct _limits
	{
		double s_epsilonSpatial;
		double s_epsilonScale;
		unsigned int s_maxNbIterAll;
		unsigned int s_maxNbIterSpatial;
		unsigned int s_maxNbIterScale;

		/**
			Copy constructor
		*/
		_limits(const _limits& source)
		: s_epsilonSpatial(source.s_epsilonSpatial), s_epsilonScale(source.s_epsilonScale),
		  s_maxNbIterAll(source.s_maxNbIterAll), s_maxNbIterSpatial(source.s_maxNbIterSpatial),
		  s_maxNbIterScale(source.s_maxNbIterScale)
		{}

		/**
			Default constructor
		*/
		_limits(){};

		/**
			Parametrized constructor
		*/
		_limits(double epsSpatial, double epsScale, unsigned int maxAll, unsigned int maxSpatial,
				unsigned int maxScale)
		: s_epsilonSpatial(epsSpatial), s_epsilonScale(epsScale), s_maxNbIterAll(maxAll),
		  s_maxNbIterSpatial(maxSpatial), s_maxNbIterScale(maxScale)
		{}

	}limits;


public:

	/**
		Default constructor
		\param model Pointer to previously acquired target model
		\param b Scale range step (usually 1.1)
		\param n Scale range radius (usually 2)
		\param bounds Convergence parameters
		\param Optional feature descriptor filter
	*/
	TrackMeanShift(TrackMSTargetBase* model, double b, int n, limits bounds,
			trackFilter* filt = nullptr, TrackMatch* match = nullptr);

	/**
		Copy constructor, useful since tracker objects are passed in std::vector (therefore are copied)
	*/
	TrackMeanShift(const TrackMeanShift& source);

public:

	/**
		Main tracking function for meanShift, parameters are returned in base data members of trackFlow.
	*/
	bool track(unsigned char* image);

	/**
		Norm function
	*/
	double norm(double loc1[2], double x, double y);

protected:

	TrackMSTargetBase*	m_model;		/*<< Pointer to meanshift target (includes derived classes): this determines the type of meanshift performed */
	TrackMSTargetBase*  m_candidate;	/*<< Pointer to current candidate */

	double m_b;							/*<< Scale range factor */
	int m_n;							/*<< Nb of scales tracked */

	limits m_limits;					/*<< Stopping parameters */

};