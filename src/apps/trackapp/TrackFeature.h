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
class   trackFeatureSIFT: public TrackFeature
{

public:

	/**
		Default constructor
	*/
	trackFeatureSIFT()
	: TrackFeature(0, 0), _desc(nullptr), _size(0)
	{}

	/**
		Constructor with descriptor array allocation and location initialization
	*/
	trackFeatureSIFT(double x, double y, unsigned int size)
	: TrackFeature(x, y), _size(size)
	{
		_desc = new vl_sift_pix[size];
	}

	/**
		Destructor, deallocates m_desc
	*/
	virtual
	~trackFeatureSIFT()
	{
		delete[] _desc;
	}

	/**
		Access histogram
	*/
	vl_sift_pix* desc();

	/**
		Comparator: Euclidean distance for base SIFT
	*/
	double compare(TrackFeature* x);

protected:
 	vl_sift_pix* _desc;	//!< Associated descriptor
 	unsigned int _size; //!< Size of descriptor array
};

inline
vl_sift_pix* trackFeatureSIFT::desc()
{
	return _desc;
}

/**
	Feature point array
*/
class   trackFeatArray
{

public:
	/**
		Destructor
	*/
	virtual
	~trackFeatArray()
	{}

	/**
		returns size
	*/
	unsigned int size() const
	{
		return _feature.size();
	}

	/**
		returns array of features
	*/
	TrackFeature* operator[](unsigned int index)
	{
		return _feature[index];
	}

	/**
		Returns the vector of features
	*/
	std::vector<TrackFeature*>& feature()
	{
		return _feature;
	}

	/**
		Show feature points in the array
	*/
	void show() ;

protected:

	std::vector<TrackFeature*> _feature; //!< Array of feature points
//s	unsigned int _size; //!< Size of the array
};

/**
	Base filter class
*/
class  trackFilter
{
public:

	/**
		Default constructor
	*/
	trackFilter()
	{}

	virtual
	~trackFilter()
	{}

public:

	/**
		To remove keypoints that are not in the model shape.
	*/
	virtual
	void prune(unsigned char* validity, unsigned int width) = 0;

	/**
		Computes feature descriptors
	*/
	virtual
	void compute(unsigned char* image, unsigned int width, unsigned int height,
			bool isModel = false, unsigned char* bw = nullptr) = 0;

	/**
		Additional feature descriptor processing
	*/
	virtual
	void process() = 0;

	/**
		Return feature array
	*/
	trackFeatArray* feature();

protected:

	trackFeatArray m_pFeature;		//!< Pointer to the feature array
};

inline
trackFeatArray* trackFilter::feature()
{
	return &m_pFeature;
}


/**
	SIFT Filter wrapper
*/
class  trackSIFTFilter: public trackFilter
{
public:

	/**
		Constructor
	*/
	trackSIFTFilter(int noctaves, int nlevels, int o_min);

	/**
		Copy constructor
	*/
	trackSIFTFilter(const trackSIFTFilter& source)
	: _filter(nullptr), _kp(nullptr), _noctaves(source._noctaves), _nlevels(source._nlevels),
	  _nomin(source._nomin)
	{}

	/**
		Destructor
	*/
	virtual
	~trackSIFTFilter();

	/**
		SIFT transform on a particular image
	*/
	void compute(unsigned char* image, unsigned int width, unsigned int height,
			bool isModel = 0, unsigned char* bw = nullptr);
	
	/**
		Process keypoints
	*/
	void process();

	/**
		To remove keypoints that are not in the model shape.
	*/
	void prune(unsigned char* validity, unsigned int width);

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
int trackSIFTFilter::octaves() const
{
	return _noctaves;
}

inline
int trackSIFTFilter::levels() const
{
	return _nlevels;
}

inline
int trackSIFTFilter::omin() const
{
	return _nomin;
}
