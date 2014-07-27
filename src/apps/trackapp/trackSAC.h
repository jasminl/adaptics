#pragma once

#ifdef TRACK_EXPORTS
#define TRACK_API __declspec(dllexport)
#else
#define TRACK_API __declspec(dllimport)
#endif

class trackMatch;

/**
	Base class for random consensus matching algorithms
*/
class  TRACK_API trackSAC
{
public: //typedefs

	/**
		Structure that defines rotation, scaling and translation transforms
	*/
	typedef struct _RST
	{
		double s_rot;		/*<< Rotation angle */
		double s_scale;		/*<< Scaling        */
		double s_t1;		/*<< Horizontal translation */
		double s_t2;		/*<< Vertical translation */

		/**
			Default constructor
		*/
		_RST(){s_rot = s_scale = s_t1 = s_t2 = 0;}

		/**
			Constructor
		*/
		_RST(double rot, double scale, double t1, double t2):s_rot(rot),s_scale(scale),s_t1(t1),s_t2(t2)
		{}

	}rst;

public:	//Member functions

	/**
		Estimates parameters for rotation, scale and translation
	*/
	virtual void estimateRST(rst& out){}	//TODO: maybe this shouldn't be derived...
	
	/**
		Estimates number of iterations to be run
	*/
	double estimateIter(double probQ);

	/**
		Estimates inlier MSS probability
	*/
	double estimateInlierMSSProb(unsigned int N, unsigned int NInlier, unsigned int k);

	/**
		Get minimal sample set 
	*/
	void mss(rst& out){};//TODO: complete this

public:

	/**
		Default constructor
	*/
	trackSAC(void)
	{
		m_match = NULL;
		m_alarmRate = 0.000001;		
		m_minIter = 100;			
		m_maxIter = 100000;
		m_adaptIter = 0;
	}

	/**
		Constructor with matched point object
	*/
	trackSAC(trackMatch* pt, double alarmRate=0.000001,
		unsigned int minIter = 100,
		unsigned int maxIter = 10000):m_match(pt),m_alarmRate(alarmRate),m_minIter(minIter),m_maxIter(maxIter)
	{
		m_adaptIter = 0;
	}

	/**
		Destructor
	*/
	virtual ~trackSAC(void){}
	
protected:	//Data members
	
	trackMatch*	m_match;			/*<<Matched points */

	double m_alarmRate;				/*<<Threshold to determine number of iterations */
	unsigned int m_maxIter;			/*<<Maximum number of iterations to perform SAC */
	unsigned int m_minIter;			/*<<Minimum number of iterations to perform SAC */

	unsigned int m_adaptIter;		/*<<Number of iterations to perform based on calculations */
};

class  TRACK_API trackRANSAC: public trackSAC
{
public:	//Member functions

	void estimateRST(rst& out);

public: //Constructors/destructor
	
	/**
		Default constructor
	*/
	trackRANSAC(){}

	/**
		Constructor with matched points
	*/
	trackRANSAC(trackMatch* pt, double alarmRate):trackSAC(pt,alarmRate)
	{}
};