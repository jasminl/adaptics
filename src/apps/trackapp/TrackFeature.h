#pragma once

#include <vector>
extern "C" {
#include "vl/sift.h"
}

/**
	Base class for feature point descriptors
*/
class  TrackFeature
{

public:
	/**
		Constructor with location
	*/
	TrackFeature(double x, double y)
	: _x(x), _y(y)
	{}

	virtual
	~TrackFeature(void)
	{}

	/**
		Function to compare two feature points (useful during matching)
	*/
	virtual
	double compare(TrackFeature* x) = 0;

	/**
		Returns spatial location
	*/
	std::pair<double,double> spatial() const;

	/**
		To display feature point
	*/
	void show() const;

protected:

	double _x;				//!< coordinate
	double _y;				//!< coordinate
};


inline
std::pair<double,double> TrackFeature::spatial() const
{
	return std::make_pair(_x, _y);
}


/**
	SIFT feature descriptor: this may be useless
*/
class TrackFeatureSIFT: public TrackFeature
{

public:

	/**
		Default constructor
	*/
	TrackFeatureSIFT()
	: TrackFeature(0, 0)
	{}

	/**
		Constructor with descriptor array allocation and location initialization
	*/
	TrackFeatureSIFT(double x, double y, unsigned int size)
	: TrackFeature(x, y)
	{
		_desc.resize(size);
	}

	/**
		Destructor
	 */
	virtual
	~TrackFeatureSIFT()
	{}

	/**
		Access histogram
	*/
	vl_sift_pix* desc();

	/**
		Comparator: Euclidean distance for base SIFT
	*/
	double compare(TrackFeature* x);

protected:
	std::vector<vl_sift_pix> _desc;
};

inline
vl_sift_pix* TrackFeatureSIFT::desc()
{
	return &_desc[0];
}

/**
	Feature point array
*/
class TrackFeatArray
{

public:

	/**
		Destructor
	*/
	virtual
	~TrackFeatArray()
	{}

	/**
		returns size
	*/
	unsigned int size() const;

	/**
		returns array of features
	*/
	TrackFeature* operator[](unsigned int index);

	/**
		Returns the vector of features
	*/
	std::vector<TrackFeature*>& feature();

	/**
		Show feature points in the array
	*/
	void show() ;

protected:

	std::vector<TrackFeature*> _feature; //!< Array of feature points
};

inline
unsigned int TrackFeatArray::size() const
{
	return _feature.size();
}

inline
TrackFeature* TrackFeatArray::operator[](unsigned int index)
{
	return _feature[index];
}

inline
std::vector<TrackFeature*>& TrackFeatArray::feature()
{
	return _feature;
}

/**
	Base filter class
*/
class  TrackFilter
{
public:

	virtual
	~TrackFilter()
	{}

public:

	/**
		To remove keypoints that are not in the model shape.
	*/
	virtual
	void prune(std::vector<unsigned char>& validity, unsigned int width) = 0;

	/**
		Computes feature descriptors
	*/
	virtual
	void compute(unsigned char* image, unsigned int width, unsigned int height,
			bool isModel, std::vector<unsigned char>& bw) = 0;

	/**
		Additional feature descriptor processing
	*/
	virtual
	void process() = 0;

	/**
		Return feature array
	*/
	TrackFeatArray* feature();

protected:

	TrackFeatArray _p_feat;		//!< Pointer to the feature array
};

inline
TrackFeatArray* TrackFilter::feature()
{
	return &_p_feat;
}


/**
	SIFT Filter wrapper
*/
class  TrackSIFTFilter: public TrackFilter
{
public:

	/**
		Constructor
	*/
	TrackSIFTFilter(int noctaves, int nlevels, int o_min);

	/**
		Copy constructor
	*/
	TrackSIFTFilter(const TrackSIFTFilter& source)
	: _filter(nullptr), _kp(nullptr), _noctaves(source._noctaves), _nlevels(source._nlevels),
	  _nomin(source._nomin)
	{}

	/**
		Destructor
	*/
	virtual
	~TrackSIFTFilter();

	/**
		SIFT transform on a particular image
	*/
	void compute(unsigned char* image, unsigned int width, unsigned int height,
			bool isModel, std::vector<unsigned char>& bw);
	
	/**
		Process keypoints
	*/
	void process();

	/**
		To remove keypoints that are not in the model shape.
	*/
	void prune(std::vector<unsigned char>& validity, unsigned int width);

	int octaves() const;

	int levels() const;

	int omin() const;

protected:

	_VlSiftFilt* _filter;
	const VlSiftKeypoint* _kp;
	int _noctaves;
	int _nlevels;
	int _nomin;
};

inline
int TrackSIFTFilter::octaves() const
{
	return _noctaves;
}

inline
int TrackSIFTFilter::levels() const
{
	return _nlevels;
}

inline
int TrackSIFTFilter::omin() const
{
	return _nomin;
}
