#include "TrackMatch.h"
#include <iostream>

using namespace std;

void trackMatchTri::match(TrackFeatArray* model, TrackFeatArray* target)
{
	int size = (int)model->size();	//Get number of points

	_x1.clear();						//Empty matched arrays
	_x2.clear();

	if (size == 0)
	{
		cout<<"No features to match in model\n";
		return;				//Means there are no features to match
	}

	for(int i = 0; i < size; i++)
	{
		auto res = firstSecond((*model)[i], target);		//Get first and second neighbors
		if(res.first.second / res.second.second < _t)
		{
			//This is a valid pair: add to list
			_x1.push_back((*model)[i]->spatial());						//Add model point
			_x2.push_back((*model)[res.first.first]->spatial());		//Add target point
		}
	}
}

pair< pair<unsigned int, double>, pair<unsigned int, double> > trackMatchTri::firstSecond(TrackFeature* modelPt, TrackFeatArray* target)
{
	double delta1 = 1000;		//Start with an artificially high value
	double delta2 = 2*delta1;
	double d = 0;

	unsigned int pt1 = 0, pt2 = 0;

	for(int i = 0;i < (int)target->size(); i++)
	{
		if((d = modelPt->compare((*target)[i])) < delta1)
		{
			//This is so far the closest neighbor
			pt1 = (unsigned)i;			//keep ID of point
			delta1 = d;			//Update delta1 using new minimal distance
		}
		else if(d < delta2)
		{
			//This is farther than first neighbor, but closer than second
			pt2 = (unsigned)i;
			delta2 = d;
		}
	}

	//Return both point IDs and distances
	return make_pair(make_pair(pt1,delta1), make_pair(pt2, delta2));
}
