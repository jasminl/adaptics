#include "trackSAC.h"

#include <cmath>

double trackSAC::estimateIter(double probQ)
{
	return ceil(log(m_alarmRate)/log(1 - probQ));
}

double trackSAC::estimateInlierMSSProb(unsigned int N, unsigned int NInlier, unsigned int k)
{
	double q=1;

	for(unsigned int i=0;i<k;i++)
	{
		q *= (NInlier-i)/(N-i);	
	}

	return q;
}

void trackRANSAC::estimateRST(rst& out)
{
	rst candidate;		//Declare temporary buffer
	
	unsigned int i=0;	//Iteration number

	while(i<m_minIter || (i<m_adaptIter && i<m_maxIter))
	{
		//Get a minimum sampling set

		//Get consensus set

		//TODO: Normalize points?

		i++;		//Add one iteration
	}
}
