#pragma once

#include <vector>
extern "C" {
#include "vl/sift.h"
}


using namespace std;

/**
	Base class for feature point descriptors
*/
class  trackFeature
{

public:	//Member functions

	/**
		Function to compare two feature points (useful during matching)
	*/
	virtual double compare(trackFeature* x)=0;

	/**
		Returns spatial location
	*/
	pair<double,double> spatial() const
	{
		return pair<double,double> (m_x,m_y);
	}

	/**
		To display feature point
	*/
	void show() const;

public:
	trackFeature(void){}
	
	/**
		Constructor with location
	*/
	trackFeature(double x,double y):m_x(x),m_y(y)
	{}

	virtual ~trackFeature(void){}

protected:	//Data member

	double m_x;				/*<<x coordinate */
	double m_y;				/*<<y coordinate */
};

/**
	SIFT feature descriptor: this may be useless
*/
class   trackFeatureSIFT: public trackFeature
{
protected:	//Data members

	unsigned int	m_size;			/*<<Size of descriptor array */
 	vl_sift_pix*    m_desc;			/*<<Associated descriptor */

public:	//Constructors /destructor

	/**
		Default constructor
	*/
	trackFeatureSIFT()
	{
		m_desc = NULL;
		m_size = 0;
	}

	/**
		Constructor with descriptor array allocation and location initialization
	*/
	trackFeatureSIFT(double x, double y,unsigned int size ):trackFeature(x,y),m_size(size)
	{
		m_desc = new vl_sift_pix[size];
	}

	/**
		Destructor, deallocates m_desc
	*/
	virtual ~trackFeatureSIFT()
	{
		delete[] m_desc;
	}

public:	//Member functions

	/**
		Access histogram
	*/
	vl_sift_pix* desc()
	{
		return m_desc;
	}

	/**
		Comparator: Euclidean distance for base SIFT
	*/
	double compare(trackFeature* x)
	{
		trackFeatureSIFT* px = dynamic_cast<trackFeatureSIFT*>(x);

		double d=0;

		for(unsigned int i=0;i<m_size;i++)
		{
			d += (m_desc[i] - px->desc()[i]) * (m_desc[i] - px->desc()[i]);
		}

		return d;
	}

};


/**
	Feature point array
*/
class   trackFeatArray
{
public:	//Member functions

	/**
		returns size
	*/
	unsigned int size() const
	{
		return m_feature.size();
	}

	/**
		returns array of features
	*/
	trackFeature* operator[](unsigned int index)
	{
		return m_feature[index];
	}

	/**
		Returns the vector of features
	*/
	vector<trackFeature*>& feature()
	{
		return m_feature;
	}

	/**
		Show feature points in the array
	*/
	void show() ;

public:	//Constructors/ destructor

	/**
		Default constructor
	*/
	trackFeatArray(){}

	/**
		Destructor
	*/
	virtual ~trackFeatArray(){}

protected:	//Data members

	vector<trackFeature*> m_feature;	/*<<Array of feature points */

	unsigned int m_size;		/*<<Size of the array */

};


/**
	Base filter class
*/
class  trackFilter
{
public:		//Member functions

	/**
		To remove keypoints that are not in the model shape.
	*/
	virtual void prune(unsigned char* validity,unsigned int width)=0;

	/**
		Computes feature descriptors
	*/
	virtual void compute(unsigned char* image, unsigned int width, unsigned int height,bool isModel=false, unsigned char* bw=NULL)=0;

	/**
		Additional feature descriptor processing
	*/
	virtual void process()=0;

	/**
		Return feature array
	*/
	trackFeatArray* feature()
	{
		return &m_pFeature;
	}

protected:	//Data member

	trackFeatArray m_pFeature;		/*<< Pointer to the feature array */

public:	//Constructors / destructor
	
	/**
		Default constructor
	*/
	trackFilter()
	{
	}


};

/**
	SIFT Filter wrapper
*/
class  trackSIFTFilter: public trackFilter
{
public:	//Member functions

	/**
		SIFT transform on a particular image
	*/
	void compute(unsigned char* image, unsigned int width, unsigned int height, bool isModel=0,unsigned char* bw=NULL);
	
	/**
		Process keypoints
	*/
	void process();

	/**
		To remove keypoints that are not in the model shape.
	*/
	void prune(unsigned char* validity,unsigned int width);

	int octaves() const
	{
		return m_noctaves;
	}

	int levels() const
	{
		return m_nlevels;
	}

	int omin() const
	{
		return m_nomin;
	}

public:

	/**
		Constructors
	*/
	trackSIFTFilter(int noctaves, int nlevels, int o_min);

	/**
		Copy constructor
	*/
	trackSIFTFilter(const trackSIFTFilter& source):m_noctaves(source.m_noctaves),m_nlevels(source.m_nlevels),m_nomin(source.m_nomin)
	{
		m_kp = NULL;
	}

	/**
		Destructor
	*/
	~trackSIFTFilter();

protected:

	_VlSiftFilt* m_filter;
	const VlSiftKeypoint* m_kp;

	int m_noctaves;
	int m_nlevels;
	int m_nomin;

	

};
