#pragma once

#include "TrackMSTargetColor.h"

/**
	RGB color mean-shift target
*/
class TrackMSTargetRGB: public TrackMSTargetColor
{
public:

	/**
	 * Constructor
	 * \param x x-coordinate of target ellipse
	 * \param y y-coordinate of target ellipse
	 * \param hx Width of target ellipse
	 * \param hy Height of target ellipse
	 * \param s Scales to track
	 * \param image Input frame
	 * \param image_width Width of the input frame
	 * \param image_height Height of the input frame
	 * \param back_size Size of background for background reweighting (if -1, no background reweighting is used)
	 * \param rNbBin Number of bins in red dimension
	 * \param gNbBin Number of bins in green dimension
	 * \param bNbBin Number of bins in blue dimension
	 */
	TrackMSTargetRGB(double x, double y, double hx, double hy,
			std::vector<double> s, unsigned char* image, int image_width, int image_height,
			double back_size = -1, unsigned int rNbBin = 16, unsigned int gNbBin = 16, unsigned int bNbBin = 16);

	/**
		Bin assignment function
	*/
	unsigned int quantize(unsigned char r,unsigned  char g,unsigned  char b);

	/**
		Create target candidate from previous candidate
	*/
	void candidate(TrackMSTargetBase*& dest, unsigned char* image, double x, double y,
			std::vector<double> scale);

	/**
		Updates target candidate without reallocating it
	*/
	void update_candidate(TrackMSTargetBase* dest, unsigned char* image, double x, double y,
			std::vector<double> scale);

	/**
		Destroys allocated candidate
	*/
	void kill_candidate(TrackMSTargetBase*& dest);
};
