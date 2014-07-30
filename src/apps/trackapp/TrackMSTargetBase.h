#pragma once

#include <vector>
#include "fastmap.h"

/**
	\brief Class for representing a mean-shift target
	Whatever feature is used, a target should have a center, scale and kernel
*/
class TrackMSTargetBase
{
public:

	/**
		Default constructor
	*/
	TrackMSTargetBase();

	/**
		Constructor used to initialize target
	*/
	TrackMSTargetBase(double x, double y, double hx, double hy,
			std::vector<double> s, unsigned char* image, unsigned int imageWidth, unsigned int imageHeight);

	/**
	 * Destructor, does nothing.
	 */
	virtual
	~TrackMSTargetBase()
	{}

	/**
		Allocates arrays
	*/
	void allocate(std::vector<double>& s);

	/**
		Create target candidate from previous candidate
	*/
	virtual
	void candidate(TrackMSTargetBase*& dest, unsigned char* image, double x, double y,
			std::vector<double> scale) = 0;

	/**
		Updates target candidate
	*/
	virtual
	void update_candidate(TrackMSTargetBase* dest, unsigned char* image, double x, double y,
			std::vector<double> scale) = 0;
	
	/**
		Destroy candidate
	*/
	virtual
	void kill_candidate(TrackMSTargetBase*& dest) = 0;

	/**
		Return current x coordinate
	*/
	double current_x() const;

	/**
		Return current y coordinate
	*/
	double current_y() const;

	/**
		Return current scale
	*/
	double current_scale() const;

	/**
		To set x, y and scale
	*/
	void set(double x, double y, const std::vector<double> s);

	/**
		Return current scale with allowed modification
	*/
	double& c_scale();

	/**
		Generate vector of all tracked scales given b, n and current center scale
	*/
	std::vector<double> expand_scale(double b, int n);

	/**
		Epanichnikov kernel
	*/
	double epanichnikov(double distance, double cd, double d);
	
	/**
		Epanichnikov kernel in scale-space version
	*/
	double hs(double scale, double n);

	/**
		Kernel function for Gaussian pyramid
	*/
	double kp(double distance, double sigma);

	/**
		Kernel function for scale mean shift	(Eq.17)
	*/
	double hxq(double distance, double sigma);

	/**
		Distance method
	*/
	virtual double bhattacharyya_distance(TrackMSTargetBase* a, double scale);

	/**
		Obtain bounding rectangle start and end points
	*/
	void bounding_rect(double verticalAxis, double horizontalAxis);
	
	/**
		Setup for scale manipulations in both spatial and scale mean shift std::vectors
	*/
	std::vector<double> setup_scale(int& n);

	/**
		Spatial meanshift std::vector
	*/
	std::pair<double, double> spatial_meanshift();

	/**
		Scale meanshift std::vector
	*/
	double scale_meanshift();

	/**
		Weight function
	*/
	virtual void weight(TrackMSTargetBase* model)=0;

	/**
		Get the scale range
	*/
	std::vector<double>& scale_range();

	/**
		Get image widht and height
	*/
	std::pair<unsigned int, unsigned int> wh() const;

protected:
	double _x;				//!< x-coordinate of the center
	double _y;				//!< y-coordinate of the center

	double _hx;				//!< width of the tracking ellipse
	double _hy;				//!< height of the tracking ellipse

	std::vector<double> _s;	//!< scale of the tracking ellipse

	unsigned char* _image;	//!< Pointer to image from which to get representation
	unsigned int _width;	//!< Image width
	unsigned int _height;   //!< Image height

	unsigned int _vertical1;	//!< Rectangle upper pixel
	unsigned int _vertical2;	//!< Rectangle lower pixel
	unsigned int _horizontal1;	//!< Rectangle left pixel
	unsigned int _horizontal2;	//!< Rectangle right pixel

	fastmap<double,double> _w2;								//!< Weight
	fastmap<double,unsigned int> _bin2;						//!< Bin assignments
	fastmap<double, std::pair<double,double> > _cenPix2;	//!< Centered pixels
	fastmap<double, double> _dist2;							//!< Pixel distances
};

inline
double TrackMSTargetBase::current_x() const
{
	return _x;
}

inline
double TrackMSTargetBase::current_y() const
{
	return _y;
}
