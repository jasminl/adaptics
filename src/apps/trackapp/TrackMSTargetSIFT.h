#pragma once

/**
	Class for representing a mean-shift SIFT target
*/
class TrackMSTargetSIFT: public TrackMSTargetBase
{
public:

	/**
		Create target candidate from previous candidate
	*/
	virtual
	void candidate(TrackMSTargetBase*& dest, unsigned char* image, double x, double y,
			std::vector<double> scale)
	{
		dest = nullptr;
	}

	/**
		Updates target candidate
	*/
	void update_candidate(TrackMSTargetBase* dest, unsigned char* image, double x, double y,
			std::vector<double> scale)
	{

	}

	/**
		Destroy candidate
	*/
	void kill_candidate(TrackMSTargetBase*& dest)
	{
		delete dynamic_cast<TrackMSTargetSIFT*>(dest);
		dest = nullptr;
	}
};
