#pragma once

#include "trackFlow.h"

/**
	Main class for particle filter tracking
*/
class  trackParticleFilter: public trackFlow
{
public:

	bool track(unsigned char* image)
	{
		return true;
	}

public:	//Constructors /destructors
	trackParticleFilter(trackFilter* filt, trackMatch* match):trackFlow(filt,match)			//Optional feature descriptor filter)
	{}

};

