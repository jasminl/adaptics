#pragma once

#include "TrackFlow.h"

/**
	Main class for particle filter tracking
*/
class  TrackParticleFilter: public TrackFlow
{

public:

	/**
	 * Base constructor
	 */
	TrackParticleFilter(trackFilter* filt, TrackMatch* match)
	: TrackFlow(filt, match)
	{}

	/**
	 * Destructor
	 */
	virtual
	~TrackParticleFilter()
	{}

	/**
	 * Returns true for now
	 * \note This function is just a placeholder
	 */
	bool track(unsigned char* image);
};

inline
bool TrackParticleFilter::track(unsigned char* image)
{
	return true;
}
