#include "TrackMeanShift.h"
#include <cmath>

inline
double TrackMeanShift::norm(double loc1[2], double x, double y)
{
	return sqrt((loc1[0] - x)*(loc1[0] - x) + (loc1[1] - y)*(loc1[1] - y));
}

TrackMeanShift::TrackMeanShift(TrackMSTargetBase* model, double b, int n, limits bounds,
		TrackFilter* filt, TrackMatch* match)
: TrackFlow(filt, match), _model(model), _candidate(nullptr), _scale_range(b), _nb_scales(n),
  _limits(bounds)
{
	//Obtain current scale-space locations (for scale, always middle scale in vector
	_cur_scale = model->current_scale();
	_cur_x = model->current_x();
	_cur_y = model->current_y();
	_cur_expanded_scale = _model->expand_scale(_scale_range, _nb_scales);
}


TrackMeanShift::TrackMeanShift(const TrackMeanShift& source):
		TrackFlow(source._feat_filter,source._match),_model(source._model), _candidate(nullptr),
		_scale_range(source._scale_range), _nb_scales(source._nb_scales), _limits(source._limits)
{}

bool TrackMeanShift::track(unsigned char* image)
{
	double s1 = 1000;            //Initial conditions for scale test in while loop
	double y1[2] = {1000,1000};  //"  "   "   "   "   "   location test " "   "
	int nbIterAll=0;             //Number of iterations of entire algorithm
	bool badLoc   = false;       //Indicator to see if location shift has converged

	//General tracking
	while(((_cur_scale - s1) >= _limits.s_epsilonScale || norm(y1, _cur_x, _cur_y) >= _limits.s_epsilonSpatial)
			&& (nbIterAll < _limits.s_maxNbIterAll))
	{

		//Track in space
		_model->expand_scale(_scale_range, _nb_scales);

		_model->candidate(_candidate,image,_cur_x,_cur_y,_cur_expanded_scale);		//Instantiation of candidate representation

		double rho0 = _candidate->bhattacharyya_distance(_model,_cur_scale);	//Get initial rho
		double rho1;

		bool shiftLeg = false, firstIteration = true;

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
				_model->candidate(_candidate,image,_cur_x,_cur_y,_cur_expanded_scale);
			}

			_candidate->weight(_model);

			auto shift = _candidate->spatial_meanshift();
			y1[0] += shift.first;
			y1[1] += shift.second;

			_model->update_candidate(_candidate,image,_cur_x,_cur_y,_cur_expanded_scale);
			rho1 = _candidate->bhattacharyya_distance(_model,_cur_scale);

			int meanShiftIter = 1;
			while(rho1<rho0)
			{
				y1[0] = 0.5*(_cur_x + y1[0]);
				y1[1] = 0.5*(_cur_y + y1[1]);

				//Calculate p1 and correlate
				_model->update_candidate(_candidate,image,_cur_x,_cur_y,_cur_expanded_scale);

				rho1 = _candidate->bhattacharyya_distance(_model,_cur_scale);

				meanShiftIter++;

				if(meanShiftIter > _limits.s_maxNbIterSpatial)
				{
					badLoc = true;
					break;
				}
			}

			if(norm(y1,_cur_x,_cur_y)<_limits.s_epsilonSpatial)
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
			_model->update_candidate(_candidate,image,_cur_x,_cur_y,_cur_expanded_scale);

			if(firstIteration==true)
			{
				//If it's the first iteration, set to initial position
				firstIteration = false;
				s1 = _cur_scale;
			}

			_candidate->weight(_model);
			double ds = _candidate->scale_meanshift();

			s1 = s1 * pow(_scale_range,ds);
			_candidate->c_scale() = s1;										//Rest scale to new location
			_cur_expanded_scale = _candidate->expand_scale(_scale_range,_nb_scales);		//Adjust searched range

			_model->update_candidate(_candidate,image,_cur_x,_cur_y,_cur_expanded_scale);
			double rho2 = _candidate->bhattacharyya_distance(_model,_candidate->current_scale());

			int meanShiftIter = 1;

			while(rho2<rho1)
			{
				s1 = 0.5 * (_cur_scale + s1);
				_candidate->c_scale() = s1;									//Rest scale to new location
				_cur_expanded_scale = _candidate->expand_scale(_scale_range,_nb_scales); //Adjust searched range

				//Calculate p1 and correlate
				_model->update_candidate(_candidate,image,_cur_x,_cur_y,_cur_expanded_scale);
				rho2 = _candidate->bhattacharyya_distance(_model,_candidate->current_scale());

				meanShiftIter++;

				if(meanShiftIter > _limits.s_maxNbIterScale)
					break;
			}

			if (sqrt((_cur_scale-s1)*(_cur_scale-s1))<_limits.s_epsilonScale)
				shiftLeg = true;


			_cur_scale = s1;	//Always set sigma0 to s1 for next iteration of entire algorithm
			_cur_expanded_scale = _candidate->expand_scale(_scale_range,_nb_scales);
		}


		nbIterAll++;
	}

	_model->kill_candidate(_candidate);

	//Transfer parameters to coord structure
	_coord = coord(_cur_x,_cur_y,_cur_scale,_cur_expanded_scale);

	//Return whether target is still good (true) or has been lost (false), note we accept lost scales
	if(_limits.s_maxNbIterAll == nbIterAll || badLoc == true)
		return false;		//Lost target
	return true;		//Valid target
}
