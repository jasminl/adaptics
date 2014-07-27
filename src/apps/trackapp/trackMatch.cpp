#include "trackMatch.h"
#include "trackFeature.h"

#include <iostream>

void trackMatchTri::match(trackFeatArray* model, trackFeatArray* target)
{
	unsigned int size = model->size();	//Get number of points

	m_x1.clear();						//Empty matched arrays
	m_x2.clear();

	if (size==0)
	{
		std::cout<<"No features to match in model\n";	
		return;				//Means there are no features to match
	}

	pair< pair<unsigned int, double>, pair<unsigned int, double> > res;

	for(unsigned int i=0;i<size;i++)
	{
		res = firstSecond((*model)[i],target);		//Get first and second neighbors
		if(res.first.second / res.second.second < m_t)
		{
			//This is a valid pair: add to list
			m_x1.push_back((*model)[i]->spatial());						//Add model point
			m_x2.push_back((*model)[res.first.first]->spatial());		//Add target point
		}
	}
}

pair< pair<unsigned int, double>, pair<unsigned int, double> > trackMatchTri::firstSecond(trackFeature* modelPt, trackFeatArray* target)
{
	unsigned int size = target->size();

	double delta1 = 1000;		//Start with an artificially high value
	double delta2 = 2*delta1;
	double d=0;

	unsigned int pt1=0, pt2=0;

	for(unsigned int i=0;i<size;i++)
	{
		if((d = modelPt->compare((*target)[i])) < delta1)
		{
			//This is so far the closest neighbor
			pt1 = i;			//keep ID of point
			delta1 = d;			//Update delta1 using new minimal distance
		}
		else if(d < delta2)
		{
			//This is farther than first neighbor, but closer than second
			pt2 = i;
			delta2 = d;
		}

	}

	//Return both point IDs and distances
	return pair< pair<unsigned int, double>, pair<unsigned int, double> >(pair<unsigned int, double>(pt1,delta1),pair<unsigned int, double>(pt2,delta2));
}
