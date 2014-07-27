// The following ifdef block is the standard way of creating macros which make exporting 
// from a DLL simpler. All files within this DLL are compiled with the TRACK_EXPORTS
// symbol defined on the command line. this symbol should not be defined on any project
// that uses this DLL. This way any other project whose source files include this file see 
// TRACK_API functions as being imported from a DLL, whereas this DLL sees symbols
// defined with this macro as being exported.

#include <vector>
#include <cmath>
#include <map>
#include <deque>
#include <iostream>

#include "trackSAC.h"

class trackFeature;			
class trackFilter;

using namespace std;

#ifdef TRACK_EXPORTS
#define TRACK_API __declspec(dllexport)
#else
#define TRACK_API __declspec(dllimport)
#endif

#define PI 3.1416

/**
	Map where key is not const
*/
template<typename T, typename Q> class fastMap
{
protected:

	vector<T> m_key;
	vector< Q* > m_data;

	vector<unsigned int> m_size;			/*<<Number of entries in the array */

	vector<unsigned int> m_maxSize;			/*<<Actual space occupied */


public:

	/**
		Constructor
	*/
	fastMap<T,Q>(vector<T> key):m_key(key){}

	/**
		Default constructor
	*/
	fastMap<T,Q>(){}
	
	~fastMap()
	{
		clear();
	}

	/**
		Add space for one element
	*/
	void add(T key, unsigned int size)
	{
		m_key.push_back(key);
		m_data.push_back(new Q[size]);
		m_size.push_back(size);

		m_maxSize.push_back(size);
	}

	/**
		Operator[]
	*/
	pair<Q*,unsigned int> operator[](T index)
	{
		vector<Q*>::iterator q = m_data.begin();
		vector<unsigned int>::iterator r = m_size.begin();

		for(vector<T>::iterator p=m_key.begin(); p != m_key.end(); p++, q++, r++)
		{
			if(*p == index)
			{
				return pair<Q*, unsigned int>(*q,*r);
			}

		}
		return pair<Q*,unsigned int>(NULL,0);	//Index not found
	}

	/**
		Operator () accesses by address
	*/
	pair<Q*, unsigned int > operator()(unsigned int index)
	{
		return pair<Q*, unsigned int >(m_data[index],m_size[index]);
	}

	/**
		Gives max size for a given scale
	*/
	unsigned& maxSize(T index)
	{
		vector<unsigned int>::iterator q = m_maxSize.begin();
		for(vector<T>::iterator p = m_key.begin() ; p != m_key.end(); p++,q++)
		{
			if(index == *p)
			{
				return *q;
			}
		}
	}

	/**
		Reallocates space and/or update index
	*/
	Q* realloc(T index, unsigned int minSize)
	{
		vector<Q*>::iterator p = m_data.begin();
		for(vector<T>::iterator q = m_key.begin() ; q != m_key.end() ; q++,p++)
		{
			if(index == *q)
			{
				if(maxSize(index) < minSize)
				{
					delete *p;

					//Reallocate only if is smaller
					*p = new Q[minSize];

					//Set to 0
					Q* arr = *p;
					memset(arr,0,sizeof(Q)*minSize);

					maxSize(index) = minSize;
					actualSize(index) = minSize;
				}
				
				return *p;	//Always exit function with current array, whether or not it was reassigned
			}
		}

		//Not found so add it
		m_key.push_back(index);
		m_size.push_back(minSize);
		m_maxSize.push_back(minSize);
		m_data.push_back(new Q[minSize]);

		return m_data[m_data.size()-1];

	}

	/**
		Change scale index
	*/
	void changeIndex(unsigned int index,T value)
	{
		if(index>=m_key.size())
		{
			m_key.push_back(value);

			m_size.push_back(0);
			m_maxSize.push_back(0);
			m_data.push_back(NULL);
		}

		m_key[index] = value;
	}

	/**
		Nb of elements
	*/
	unsigned int size() const
	{
		return (unsigned int)m_key.size();
	}

	/**
		Empty
	*/
	void clear()
	{
		m_key.empty();
		m_size.empty();
		for(vector<Q*>::iterator p = m_data.begin() ; p != m_data.end() ; p++)
		{
			delete *p;
		}
	}

	pair<Q*,unsigned int> centerScale()
	{
		unsigned int index = static_cast<unsigned int>(m_key.size()/ 2);
		return pair<Q*,unsigned int> (m_data[index],m_size[index]);

	}
 
	unsigned int& actualSize(T index)
	{
		vector<unsigned int>::iterator q = m_size.begin();
		for(vector<T>::iterator p = m_key.begin() ; p != m_key.end() ; p++, q++)
		{
			if(*p == index)
			{
				return *q;
			}
		}
	}

};
/***** Object representation classes *****/
/**
	Class for representing a mean-shift target
*/
class TRACK_API trackMSTargetBase
{

public:	//Constructors and destructors

	/**
		Default constructor
	*/
	trackMSTargetBase(){};

	/**
		Constructor used to initialize target
	*/
	trackMSTargetBase(double x, 
					  double y,
					  double hx,
					  double hy,
					  vector<double> s,
					  unsigned char* image,
					  unsigned int imageWidth,
					  unsigned int imageHeight):m_x(x),m_y(y),m_hx(hx),m_hy(hy),m_s(s),m_image(image),m_width(imageWidth),m_height(imageHeight)
	{
		allocate(s);
	};

public:	//Member functions

	/**
		Allocates arrays
	*/
	void allocate(vector<double>& s)
	{
		for(unsigned int i=0;i<m_bin2.size() ; i++)
		{
			m_bin2.changeIndex(i,s[i]);
			m_w2.changeIndex(i,s[i]);
			m_cenPix2.changeIndex(i,s[i]);
			m_dist2.changeIndex(i,s[i]);
		}
	}

	/**
		Create target candidate from previous candidate
	*/
	virtual void candidate(trackMSTargetBase*& dest, unsigned char* image, double x, double y, vector<double> scale)=0;

	/**
		Updates target candidate
	*/
	virtual void updateCandidate(trackMSTargetBase* dest, unsigned char* image, double x, double y, vector<double> scale)=0;
	
	/**
		Destroy candidate
	*/
	virtual void killCandidate(trackMSTargetBase*& dest)=0;

	/**
		Return current x coordinate
	*/
	double currentX() const
	{
		return m_x;
	}

	/**
		Return current y coordinate
	*/
	double currentY() const
	{
		return m_y;
	}
	/**
		Return current scale
	*/
	double currentScale() const
	{
		unsigned int index = static_cast<unsigned int>(m_s.size()/2);
		return m_s[index];
	}

	/**
		To set x, y and scale
	*/
	void set(double x, double y, vector<double> s)
	{
		m_x = x;
		m_y = y;
		m_s = s;
	}

	/**
		Return current scale with allowed modification
	*/
	double& cScale()
	{
		unsigned int index = static_cast<unsigned int>(m_s.size()/2);
		return m_s[index];
	}

	/**
		Generate vector of all tracked scales given b, n and current center scale
	*/
	vector<double> expandScale(double b, int n)
	{
		double c = currentScale();		//Obtain current scale (if there is only one entry, returns it)

		m_s = vector<double>(2*n+1);	//Allocate space for new scale vector
		n = -n;
		for(vector<double>::iterator p= m_s.begin(); p != m_s.end(); p++,n++)
		{
			*p = c * pow(b,n);
		}

		return m_s;
	}

	/**
		Epanichnikov kernel
	*/
	double epanichnikov(double distance, double cd, double d)
	{
		if(distance<=1)
			return 0.5 / cd * (d + 2)*(1 - distance);		
		else
			return 0;
	}
	
	/**
		Epanichnikov kernel in scale-space version
	*/
	double hs(double scale, double n)
	{
		if (n==0)
			return 1;
		else
			return 1 - (scale/n)*(scale/n);
	}

	/**
		Kernel function for Gaussian pyramid
	*/
	double kp(double distance, double sigma)
	{
		double s1 = sigma / 1.6;
		double s2 = sigma * 1.6;

		double term1 = exp((- distance * distance)/(2 * s1 * s1)) / (2 * PI * pow(s1,4.0));
		double term2 = exp((- distance * distance)/(2 * s2 * s2)) / (2 * PI * pow(s2,4.0));

		return term1 - term2; 
	}

	/**
		Kernel function for scale mean shift	(Eq.17)
	*/
	double hxq(double distance, double sigma)
	{
		double fd = 2*PI*sigma*sigma / 1.6;	//Scale divided by 1.6
		double ft = 2*PI*sigma*sigma * 1.6;	//Scale times 1.6
		
		double fdex = 2*sigma*sigma / 1.6;	//Scale divided by 1.6
		double ftex = 2*sigma*sigma * 1.6;	//Scale times 1.6

		double term1 = exp(-distance*distance / fdex ) / fd;
		double term2 = exp(-distance*distance / ftex ) / ft;

		return term1 - term2;
	}

	/**
		Distance method
	*/
	virtual double bhattacharyyaDistance(trackMSTargetBase* a, double scale)
	{
		return 0;
	}

	/**
		Obtain bounding rectangle start and end points
	*/
	void boundingRect(double verticalAxis, double horizontalAxis)
	{
		m_vertical1   = static_cast<unsigned int>(max(m_y - verticalAxis,0.0));
		m_vertical2   = static_cast<unsigned int>(min(m_y + verticalAxis,m_height-1.0));
		m_horizontal1 = static_cast<unsigned int>(max(m_x - horizontalAxis,0.0));
		m_horizontal2 = static_cast<unsigned int>(min(m_x + horizontalAxis,m_width-1.0));
	}

	
	/**
		Setup for scale manipulations in both spatial and scale mean shift vectors
	*/
	vector<double> setupScale(int& n)
	{
		vector<double> ss;

		if(m_s.size() == 1)
		{
			//Mean shift is only done in space, not scale
			ss.push_back(0);
			n = 0;
		}
		else
		{
			//Mean shift is performed in scale space
			n  = static_cast<int>(m_s.size() / 2);	//Get radius of tracked scales
			ss = vector<double>(m_s.size());

			double pn = -n;
			for(vector<double>::iterator p = ss.begin(); p != ss.end(); p++)
			{
				*p = pn++;	
			}
		}

		return ss;
	}

	/**
		Spatial meanshift vector
	*/
	pair<double, double> spatialMeanShift();

	/**
		Scale meanshift vector
	*/
	double scaleMeanShift();

	/**
		Weight function
	*/
	virtual void weight(trackMSTargetBase* model)=0;

	/**
		Get the scale range
	*/
	vector<double>& scaleRange() 
	{
		return m_s;
	}

	/**
		Get image widht and height
	*/
	pair<unsigned int, unsigned int> wh() const
	{
		return pair<unsigned int,unsigned int>(m_width,m_height);
	}

protected:

	/* Whatever feature is used, a target should have a center, scale and kernel */

	double m_x;				/**< x-coordinate of the center */
	double m_y;				/**< y-coordinate of the center */

	double m_hx;			/**< width of the tracking ellipse */
	double m_hy;			/**< height of the tracking ellipse */

	vector<double> m_s;		/**< scale of the tracking ellipse */

	unsigned char*  m_image;/**< Pointer to image from which to get representation */
	unsigned int m_width;	/**< Image width  */
	unsigned int m_height;  /**< Image height */

	/*<< Bounding rectangle indices */
	unsigned int m_vertical1;
	unsigned int m_vertical2;
	unsigned int m_horizontal1;
	unsigned int m_horizontal2;

	fastMap<double,double> m_w2;							/**<Weight */
	fastMap<double,unsigned int> m_bin2;					/**<Bin assignments */
	fastMap<double, pair<double,double> > m_cenPix2;		/**<Centered pixels*/
	fastMap<double, double> m_dist2;						/**<Pixel distances */

};


/**
	Class for representing a mean-shift COLOR histogram target
*/
class TRACK_API trackMSTargetColor: public trackMSTargetBase
{
protected:	//Data members

	fastMap<double,double> m_hist;

	map<double, vector<double> > m_histogram;	/**< Color histogram, one per scale */
	vector<double> m_backHistogram;				/**< Background color histogram (optional) */

	vector<unsigned int> m_nbBinsPerDim;		/**< The number of bins per color dimension */
	vector<double> m_nbLevelsPerBin;			/**< The corresponding number of color levels in a single bin, per dimension */

	double m_backSize;							/**< Areal factor that determines the size of the background area (relative to target) to use for background weighting */

public:	//Member functions
	
	/**
		Returns histogram parameters
	*/
	pair< vector<unsigned int> , vector<double> > histParams() const
	{
		return pair< vector<unsigned int> , vector<double> >(m_nbBinsPerDim,m_nbLevelsPerBin);
	}

	/**
		Weight for target pixels
	*/
	void weight(trackMSTargetBase* model);

	/**
		Metric based on Bhattacharyya coefficient
	*/
	double bhattacharyyaDistance(trackMSTargetBase* a, double scale);

	/**
		Bin assignment function. For color models that don't have 3 components, set the third one to zero
	*/
	virtual unsigned int quantize(unsigned char q1,unsigned char q2,unsigned char q3)=0;

	/**
		Allocate histogram space
	*/
	void allocateHistograms(vector<double> scale, bool backgroundWeight);

	/**
		Histogram building function, assumes interleaved image
	*/
	void makeHistogram(vector<double>& scale);
 
	 
	/**
		Create target candidate from previous candidate
	*/
	virtual void candidate(trackMSTargetBase*& dest, unsigned char* image, double x, double y, vector<double> scale)=0;
	
	/**
		Updates target candidate
	*/
	virtual void updateCandidate(trackMSTargetBase* dest, unsigned char* image, double x, double y, vector<double> scale)=0;
	
	/**
		Destroy candidate
	*/
	virtual void killCandidate(trackMSTargetBase*& dest)=0;

	/**
		Sets histogram parameters
	*/
	void setHistogram(vector<unsigned int> nb, vector<double> sz)
	{
		m_nbBinsPerDim = nb;
		m_nbLevelsPerBin = sz;
	}

public: //Constructors/destructor

	/**
		Default constructor for color target representation
	*/
	trackMSTargetColor(double x, 
					  double y,
					  double hx,
					  double hy,
					  vector<double> s,
					  unsigned char* image,
					  unsigned int imageWidth,
					  unsigned int imageHeight,
					  double backSize):trackMSTargetBase(x,y,hx,hy,s,image,imageWidth,imageHeight),m_backSize(backSize)
	{}
};


/**
	RGB color mean-shift target
*/
class TRACK_API trackMSTargetRGB: public trackMSTargetColor
{
public:	//Member functions

	/**
		Bin assignment function
	*/
	unsigned int quantize(unsigned char r,unsigned  char g,unsigned  char b)
	{
		unsigned int red   = static_cast<unsigned int>(floor(static_cast<double>(r)/m_nbLevelsPerBin[0]));
		unsigned int green = static_cast<unsigned int>(floor(static_cast<double>(g)/m_nbLevelsPerBin[1]));
		unsigned int blue  = static_cast<unsigned int>(floor(static_cast<double>(b)/m_nbLevelsPerBin[2]));

		return red + green * m_nbBinsPerDim[0] + blue * m_nbBinsPerDim[0]*m_nbBinsPerDim[1];
	}
	
	/**
		Create target candidate from previous candidate
	*/
	void candidate(trackMSTargetBase*& dest, unsigned char* image, double x, double y, vector<double> scale)
	{
		killCandidate(dest);
		dest = new trackMSTargetRGB(x,y,m_hx,m_hy,scale,image,m_width,m_height,m_backSize,m_nbBinsPerDim[0],m_nbBinsPerDim[1],m_nbBinsPerDim[2]);
	}

	/**
		Updates target candidate without reallocating it
	*/
	void updateCandidate(trackMSTargetBase* dest, unsigned char* image, double x, double y, vector<double> scale)
	{
		trackMSTargetRGB* p = dynamic_cast<trackMSTargetRGB*>(dest);

		pair< vector<unsigned int> , vector<double> > q = p->histParams(); 

		if(m_nbBinsPerDim != q.first)				//Check that histogram doesn't change
		{
			p->setHistogram(m_nbBinsPerDim,m_nbLevelsPerBin);
		}

		p->allocateHistograms(scale,(m_backSize>0)?true:false);	//Always reallocate to keep track of scale range change
		
		p->allocate(scale);							//Allocate remaining base histograms

		p->set(x,y,scale);							//Set new space scale location
		
		p->makeHistogram(scale);					//Make histogram
	}
	
	/**
		Destroys allocated candidate
	*/
	void killCandidate(trackMSTargetBase*& dest)
	{
		delete dynamic_cast<trackMSTargetRGB*>(dest);
		dest = NULL;
	}
 
public:	//Constructors/destructor

	trackMSTargetRGB(double x,							//x-coordinate of target ellipse center
					  double y,							//y-coordinate
					  double hx,						//Width of ellipse
					  double hy,						//Height
					  vector<double> s,					//Scales to track
					  unsigned char* image,				//Input frame
					  unsigned int imageWidth,			//Width of input frame
					  unsigned int imageHeight,			//Height
					  double backSize = -1,				//Size of background for background reweighting (if -1, no background reweighting is used)
					  unsigned int rNbBin = 16,			//Number of bins in 'red' dimension
					  unsigned int gNbBin = 16,			//Number of bins in 'green' dimension
					  unsigned int bNbBin = 16);		//Number of bins in 'blue' dimension

public: //Operators

};

/**
	Class for representing a mean-shift SIFT target
*/
class TRACK_API trackMSTargetSIFT: public trackMSTargetBase
{
public:

	/**
		Create target candidate from previous candidate
	*/
	virtual void candidate(trackMSTargetBase*& dest, unsigned char* image, double x, double y, vector<double> scale)
	{
		dest = NULL;
	}

	/**
		Updates target candidate
	*/
	void updateCandidate(trackMSTargetBase* dest, unsigned char* image, double x, double y, vector<double> scale)
	{

	}
	
	/**
		Destroy candidate
	*/
	void killCandidate(trackMSTargetBase*& dest)
	{
		delete dynamic_cast<trackMSTargetSIFT*>(dest);
		dest = NULL;
	}
};

/***** Tracker classes *****/

/**
	Base class for all tracking methods
*/
class TRACK_API trackFlow
{
public: //Typedefs

	/**
		Structure used to retrieve parameters
	*/
	typedef struct _coord
	{
		double s_x;
		double s_y;
		double s_scale;
		vector<double> s_scaleRange;

		/**
			Default constructor
		*/
		_coord(){}

		/**
			Constructor
		*/
		_coord(double x, double y, double scale, vector<double> scaleRange):s_x(x),s_y(y),s_scale(scale),s_scaleRange(scaleRange)
		{}

		/**
			Display parameters
		*/
		void show() const
		{
			cout<<"(X,Y)=("<<s_x<<","<<s_y<<") Scale="<<s_scale<<" Range=(";
			
			for(unsigned int i=0;i<s_scaleRange.size(); i++) cout<<s_scaleRange[i]<<" ";
			
			cout<<")\n";
		}

	}coord;

public:	//Member functions

	virtual bool track(unsigned char* image)
	{
		return true;
	}

	coord currentCoordinates() const
	{
		return m_coord;
	}

	trackFilter* featFilter()
	{
		return m_featFilter;
	}

	trackFilter*& targetFilter()
	{
		return m_targetFilter;
	}

	trackMatch* matcher()
	{
		return m_match;
	}

public:	//Constructors/destructors
	trackFlow(trackFilter* filt, trackMatch* match):m_currentX(0),m_currentY(0),m_currentScale(0)
	{
		m_featFilter = filt;
		m_targetFilter = filt;	//TODO: why is this the same as previous?
		m_match = match;
	}

protected:	//Data members
	
	/*<< The following variables are returned by tracking methods */
	double m_currentX;							/*<< Current x-coordinates of center of target (for some methods only) */
	double m_currentY;							/*<<"	"	 y-"	"	"	"	"	"	"	"	"	"	"	"	"	"  */
	double m_currentScale;						/*<<:	:	scale	"	"	"	"	"	"	"	"	"	"	"	"	"  */
	
	trackFilter*   m_featFilter;				/*<< Feature filter for the model */
	trackFilter*   m_targetFilter;				/*<< Feature filter for a particular target */
	trackMatch*    m_match;						/*<<Matching object */

	vector<double> m_currentExpandedScale;		/*<< Current range of scales tracked */

	coord m_coord;								/*<< To return coordinates */

};


/**
	Main class for particle filter tracking
*/
class TRACK_API trackParticleFilter: public trackFlow
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


/**
	Main class for meanshift tracking
*/
class TRACK_API trackMeanShift: public trackFlow
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
		_limits(const _limits& source):s_epsilonSpatial(source.s_epsilonSpatial),
									   s_epsilonScale(source.s_epsilonScale),
									   s_maxNbIterAll(source.s_maxNbIterAll),
									   s_maxNbIterSpatial(source.s_maxNbIterSpatial),
									   s_maxNbIterScale(source.s_maxNbIterScale)
		{}

		/**
			Default constructor
		*/
		_limits(){};

		/**
			Parametrized constructor
		*/
		_limits(double epsSpatial, 
			    double epsScale, 
				unsigned int maxAll,
				unsigned int maxSpatial, 
				unsigned int maxScale):s_epsilonSpatial(epsSpatial),s_epsilonScale(epsScale),s_maxNbIterAll(maxAll),s_maxNbIterSpatial(maxSpatial),s_maxNbIterScale(maxScale)
		{}

	}limits;

public:	//Member functions
	
	/**
		Main tracking function for meanShift, parameters are returned in base data members of trackFlow.
	*/
	bool track(unsigned char* image);

	/**
		Norm function
	*/
	double norm(double loc1[2], double x, double y)
	{
		return sqrt((loc1[0] - x)*(loc1[0] - x) + (loc1[1] - y)*(loc1[1] - y));
	}

public: //Constructors/destructor
	
	/**
		Default constructor
	*/
	trackMeanShift(trackMSTargetBase* model,			//Pointer to previously acquired target model
		double b,										//Scale range step (usually 1.1)
		int n,											//Scale range radius (usually 2)
		limits bounds,									//Convergence parameters
		trackFilter* filt=NULL,							//Optional feature descriptor filter
		trackMatch* match=NULL):trackFlow(filt,match),	//Optional matching object		
					   m_model(model),
					   m_b(b),
					   m_n(n),
					   m_limits(bounds)
	{

		//Obtain current scale-space locations (for scale, always middle scale in m_s vector of the model)
		m_currentScale = model->currentScale();
		m_currentX     = model->currentX();
		m_currentY     = model->currentY();
		m_candidate    = NULL;

		m_currentExpandedScale = m_model->expandScale(m_b,m_n);
	}

	/**
		Copy constructor, useful since tracker objects are passed in vector (therefore are copied)
	*/
	trackMeanShift(const trackMeanShift& source):trackFlow(source.m_featFilter,source.m_match),m_model(source.m_model),
												 m_b(source.m_b),
												 m_n(source.m_n),
												 m_limits(source.m_limits)
	{
		m_currentScale = source.m_currentScale;
		m_currentX     = source.m_currentX;
		m_currentY     = source.m_currentY;

		m_candidate = NULL;
	}

protected:	//Data members
	
	trackMSTargetBase*	m_model;		/*<< Pointer to meanshift target (includes derived classes): this determines the type of meanshift performed */
	trackMSTargetBase*  m_candidate;	/*<< Pointer to current candidate */
	
	double m_b;							/*<< Scale range factor */
	int m_n;							/*<< Nb of scales tracked */

	limits m_limits;					/*<< Stopping parameters */

};


/***** Misc functions *****/
/**
	Function to track 1 frame
*/
TRACK_API std::vector<bool> track1Frame(pair<double,double>& location,			//Output location
						   double& scale,										//Output scale
						   unsigned char* image,								//image frame
						   vector<trackFlow*> tracker,							//trackFlow objects
						   double hx,											//Width of model
						   double hy,											//Height of model
						   unsigned int width,									//Width of image
						   unsigned int height);								//Height of image


/**
	Function to crop image to region of interest, for feature matching (ASSUMES INTERLEAVED RGB)
*/

TRACK_API unsigned char* crop(unsigned char*& buffer,							//Buffer to store cropped result
							  unsigned int width,								//Desired width of crop
							  unsigned int height,								//Desired height "	"
							  double x,											//x-coordinate of center of part to crop
							  double y,											//y-coordinate of center of "	"	"
							  unsigned char* image,								//Input image from where to crop
							  unsigned int imageWidth,							//Image width
							  unsigned int imageHeight);						//Image height