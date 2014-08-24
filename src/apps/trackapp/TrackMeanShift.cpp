#include "TrackMeanShift.h"
#include <iostream> //todo remove all cout
#include <cmath>
using namespace std;

double TrackMeanShift::norm(double loc1[2], double x, double y)
{
	return sqrt((loc1[0] - x)*(loc1[0] - x) + (loc1[1] - y)*(loc1[1] - y));
}

TrackMeanShift::TrackMeanShift(TrackMSTargetBase* model, double b, int n, limits bounds,
		trackFilter* filt, TrackMatch* match)
: TrackFlow(filt, match), m_model(model), m_b(b), m_n(n), m_limits(bounds)
{
	//Obtain current scale-space locations (for scale, always middle scale in m_s std::vector of the model)
	_cur_scale = model->current_scale();
	_cur_x     = model->current_x();
	_cur_y     = model->current_y();
	m_candidate    = nullptr;
	_cur_expanded_scale = m_model->expand_scale(m_b,m_n);
}


TrackMeanShift::TrackMeanShift(const TrackMeanShift& source):TrackFlow(source._feat_filter,source._match),m_model(source.m_model),
											 m_b(source.m_b),
											 m_n(source.m_n),
											 m_limits(source.m_limits)
{
	_cur_scale = source._cur_scale;
	_cur_x     = source._cur_x;
	_cur_y     = source._cur_y;

	m_candidate = nullptr;
}



bool TrackMeanShift::track(unsigned char* image)
{
	double s1 = 1000;            //Initial conditions for scale test in while loop
	double y1[2] = {1000,1000};  //"  "   "   "   "   "   location test " "   "
	unsigned int nbIterAll=0;    //Number of iterations of entire algorithm
	bool badLoc   = false;       //Indicator to see if location shift has converged
	bool badScale = false;       //"  "   "   "   "   " scale "   "   "   "   "   "

	//General tracking
	while(((_cur_scale-s1)>=m_limits.s_epsilonScale || norm(y1,_cur_x,_cur_y)>=m_limits.s_epsilonSpatial) && (nbIterAll<m_limits.s_maxNbIterAll))
	{

		//Track in space
		m_model->expand_scale(m_b, m_n);

		m_model->candidate(m_candidate,image,_cur_x,_cur_y,_cur_expanded_scale);		//Instantiation of candidate representation

		double rho0 = m_candidate->bhattacharyya_distance(m_model,_cur_scale);	//Get initial rho
		double rho1;

		bool shiftLeg=false, firstIteration=true;

		while(shiftLeg==0)
		{
			if(firstIteration==true)
			{
				//If it's the first iteration, set to initial position
				firstIteration = false;
				y1[0] = _cur_x;
				y1[1] = _cur_y;
			}
			else
			{
				//Calculate target candidate, without recomputing rho, since it hasn't changed
				m_model->candidate(m_candidate,image,_cur_x,_cur_y,_cur_expanded_scale);
			}

			m_candidate->weight(m_model);

			pair<double,double> shift = m_candidate->spatial_meanshift();
			y1[0] += shift.first;
			y1[1] += shift.second;

			m_model->update_candidate(m_candidate,image,_cur_x,_cur_y,_cur_expanded_scale);
			rho1 = m_candidate->bhattacharyya_distance(m_model,_cur_scale);

			unsigned int meanShiftIter=1;
			while(rho1<rho0)
			{
				y1[0] = 0.5*(_cur_x + y1[0]);
				y1[1] = 0.5*(_cur_y + y1[1]);

				//Calculate p1 and correlate
				m_model->update_candidate(m_candidate,image,_cur_x,_cur_y,_cur_expanded_scale);

				rho1 = m_candidate->bhattacharyya_distance(m_model,_cur_scale);

				meanShiftIter++;

				if(meanShiftIter > m_limits.s_maxNbIterSpatial)
				{
					badLoc = true;
					std::cout<<"Bad location\n";
					break;
				}
			}

			if(norm(y1,_cur_x,_cur_y)<m_limits.s_epsilonSpatial)
			{
				shiftLeg = true;
			}

			_cur_x = y1[0];		//Always set y0 to y1 for next iteration of entire algorithm
			_cur_y = y1[1];
		}

		//Track in scale
		shiftLeg=false;
		firstIteration=true;

		while(shiftLeg==false)
		{
			m_model->update_candidate(m_candidate,image,_cur_x,_cur_y,_cur_expanded_scale);

			if(firstIteration==true)
			{
				//If it's the first iteration, set to initial position
				firstIteration = false;
				s1 = _cur_scale;
			}

			m_candidate->weight(m_model);
			double ds = m_candidate->scale_meanshift();

			s1 = s1 * pow(m_b,ds);
			m_candidate->c_scale() = s1;										//Rest scale to new location
			_cur_expanded_scale = m_candidate->expand_scale(m_b,m_n);		//Adjust searched range

			m_model->update_candidate(m_candidate,image,_cur_x,_cur_y,_cur_expanded_scale);
			double rho2 = m_candidate->bhattacharyya_distance(m_model,m_candidate->current_scale());

			unsigned int meanShiftIter = 1;

			while(rho2<rho1)
			{
				s1 = 0.5 * (_cur_scale + s1);
				m_candidate->c_scale() = s1;									//Rest scale to new location
				_cur_expanded_scale = m_candidate->expand_scale(m_b,m_n); //Adjust searched range

				//Calculate p1 and correlate
				m_model->update_candidate(m_candidate,image,_cur_x,_cur_y,_cur_expanded_scale);
				rho2 = m_candidate->bhattacharyya_distance(m_model,m_candidate->current_scale());

				meanShiftIter++;

				if(meanShiftIter> m_limits.s_maxNbIterScale)
				{
					badScale = true;
					cout<<"Bad scale\n";
					break;
				}
			}

			if (sqrt((_cur_scale-s1)*(_cur_scale-s1))<m_limits.s_epsilonScale)
			{
				shiftLeg = true;
			}

			_cur_scale = s1;	//Always set sigma0 to s1 for next iteration of entire algorithm
			_cur_expanded_scale = m_candidate->expand_scale(m_b,m_n);
		}


		nbIterAll++;
	}

	m_model->kill_candidate(m_candidate);

	//Transfer parameters to coord structure
	m_coord = coord(_cur_x,_cur_y,_cur_scale,_cur_expanded_scale);


	//Return whether target is still good (true) or has been lost (false), note we accept lost scales
	if(m_limits.s_maxNbIterAll == nbIterAll || badLoc == true)
	{
		return false;		//Lost target
	}
	else
	{
		return true;		//Valid target
	}
}