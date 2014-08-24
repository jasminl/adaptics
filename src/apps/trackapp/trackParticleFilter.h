#pragma once

#include "TrackFlow.h"

/**
	Main class for particle filter tracking
*/
class  trackParticleFilter: public TrackFlow
{
public:

	bool track(unsigned char* image)
	{
		return true;
	}

public:	//Constructors /destructors
	trackParticleFilter(trackFilter* filt, TrackMatch* match):TrackFlow(filt,match)			//Optional feature descriptor filter)
	{}

};

