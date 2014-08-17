#include "TrackSAC.h"
#include <cmath>

void TrackRANSAC::estimate_RST(rst& out)
{
	rst candidate;		//Declare temporary buffer
	
	unsigned int i = 0;	//Iteration number

	while(i < _min_iter || (i < _adapt_iter && i < _max_iter))
	{
		//Get a minimum sampling set

		//Get consensus set

		//TODO: Normalize points?

		i++;		//Add one iteration
	}
}

double TrackSAC::estimate_iter(double probQ)
{
	return ceil(log(_alarm_rate)/log(1 - probQ));
}

double TrackSAC::estimate_inlier_MSS(unsigned int N, unsigned int NInlier, unsigned int k)
{
	double q = 1;
	for(unsigned int i = 0; i < k; i++)
		q *= (NInlier - i)/(N - i);
	return q;
}
