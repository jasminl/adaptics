#include "trackMeanShift.h"
#include <iostream> //todo remove all cout
#include <cmath>
using namespace std;

double trackMeanShift::norm(double loc1[2], double x, double y)
{
	return sqrt((loc1[0] - x)*(loc1[0] - x) + (loc1[1] - y)*(loc1[1] - y));
}


/**
	Default constructor
*/
trackMeanShift::trackMeanShift(TrackMSTargetBase* model,			//Pointer to previously acquired target model
	double b,										//Scale range step (usually 1.1)
	int n,											//Scale range radius (usually 2)
	limits bounds,									//Convergence parameters
	trackFilter* filt,							//Optional feature descriptor filter
	trackMatch* match):trackFlow(filt,match),	//Optional matching object
				   m_model(model),
				   m_b(b),
				   m_n(n),
				   m_limits(bounds)
{

	//Obtain current scale-space locations (for scale, always middle scale in m_s std::vector of the model)
	m_currentScale = model->current_scale();
	m_currentX     = model->current_x();
	m_currentY     = model->current_y();
	m_candidate    = NULL;

	m_currentExpandedScale = m_model->expand_scale(m_b,m_n);
}


trackMeanShift::trackMeanShift(const trackMeanShift& source):trackFlow(source.m_featFilter,source.m_match),m_model(source.m_model),
											 m_b(source.m_b),
											 m_n(source.m_n),
											 m_limits(source.m_limits)
{
	m_currentScale = source.m_currentScale;
	m_currentX     = source.m_currentX;
	m_currentY     = source.m_currentY;

	m_candidate = nullptr;
}



bool trackMeanShift::track(unsigned char* image)
{
	double s1 = 1000;            //Initial conditions for scale test in while loop
	double y1[2] = {1000,1000};  //"  "   "   "   "   "   location test " "   "
	unsigned int nbIterAll=0;    //Number of iterations of entire algorithm
	bool badLoc   = false;       //Indicator to see if location shift has converged
	bool badScale = false;       //"  "   "   "   "   " scale "   "   "   "   "   "

	//General tracking
	while(((m_currentScale-s1)>=m_limits.s_epsilonScale || norm(y1,m_currentX,m_currentY)>=m_limits.s_epsilonSpatial) && (nbIterAll<m_limits.s_maxNbIterAll))
	{

		//Track in space
		m_model->expand_scale(m_b, m_n);

		m_model->candidate(m_candidate,image,m_currentX,m_currentY,m_currentExpandedScale);		//Instantiation of candidate representation

		double rho0 = m_candidate->bhattacharyya_distance(m_model,m_currentScale);	//Get initial rho
		double rho1;

		bool shiftLeg=false, firstIteration=true;

		while(shiftLeg==0)
		{
			if(firstIteration==true)
			{
				//If it's the first iteration, set to initial position
				firstIteration = false;
				y1[0] = m_currentX;
				y1[1] = m_currentY;
			}
			else
			{
				//Calculate target candidate, without recomputing rho, since it hasn't changed
				m_model->candidate(m_candidate,image,m_currentX,m_currentY,m_currentExpandedScale);
			}

			m_candidate->weight(m_model);

			pair<double,double> shift = m_candidate->spatial_meanshift();
			y1[0] += shift.first;
			y1[1] += shift.second;

			m_model->update_candidate(m_candidate,image,m_currentX,m_currentY,m_currentExpandedScale);
			rho1 = m_candidate->bhattacharyya_distance(m_model,m_currentScale);

			unsigned int meanShiftIter=1;
			while(rho1<rho0)
			{
				y1[0] = 0.5*(m_currentX + y1[0]);
				y1[1] = 0.5*(m_currentY + y1[1]);

				//Calculate p1 and correlate
				m_model->update_candidate(m_candidate,image,m_currentX,m_currentY,m_currentExpandedScale);

				rho1 = m_candidate->bhattacharyya_distance(m_model,m_currentScale);

				meanShiftIter++;

				if(meanShiftIter > m_limits.s_maxNbIterSpatial)
				{
					badLoc = true;
					std::cout<<"Bad location\n";
					break;
				}
			}

			if(norm(y1,m_currentX,m_currentY)<m_limits.s_epsilonSpatial)
			{
				shiftLeg = true;
			}

			m_currentX = y1[0];		//Always set y0 to y1 for next iteration of entire algorithm
			m_currentY = y1[1];
		}

		//Track in scale
		shiftLeg=false;
		firstIteration=true;

		while(shiftLeg==false)
		{
			m_model->update_candidate(m_candidate,image,m_currentX,m_currentY,m_currentExpandedScale);

			if(firstIteration==true)
			{
				//If it's the first iteration, set to initial position
				firstIteration = false;
				s1 = m_currentScale;
			}

			m_candidate->weight(m_model);
			double ds = m_candidate->scale_meanshift();

			s1 = s1 * pow(m_b,ds);
			m_candidate->c_scale() = s1;										//Rest scale to new location
			m_currentExpandedScale = m_candidate->expand_scale(m_b,m_n);		//Adjust searched range

			m_model->update_candidate(m_candidate,image,m_currentX,m_currentY,m_currentExpandedScale);
			double rho2 = m_candidate->bhattacharyya_distance(m_model,m_candidate->current_scale());

			unsigned int meanShiftIter = 1;

			while(rho2<rho1)
			{
				s1 = 0.5 * (m_currentScale + s1);
				m_candidate->c_scale() = s1;									//Rest scale to new location
				m_currentExpandedScale = m_candidate->expand_scale(m_b,m_n); //Adjust searched range

				//Calculate p1 and correlate
				m_model->update_candidate(m_candidate,image,m_currentX,m_currentY,m_currentExpandedScale);
				rho2 = m_candidate->bhattacharyya_distance(m_model,m_candidate->current_scale());

				meanShiftIter++;

				if(meanShiftIter> m_limits.s_maxNbIterScale)
				{
					badScale = true;
					cout<<"Bad scale\n";
					break;
				}
			}

			if (sqrt((m_currentScale-s1)*(m_currentScale-s1))<m_limits.s_epsilonScale)
			{
				shiftLeg = true;
			}

			m_currentScale = s1;	//Always set sigma0 to s1 for next iteration of entire algorithm
			m_currentExpandedScale = m_candidate->expand_scale(m_b,m_n);
		}


		nbIterAll++;
	}

	m_model->kill_candidate(m_candidate);

	//Transfer parameters to coord structure
	m_coord = coord(m_currentX,m_currentY,m_currentScale,m_currentExpandedScale);


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
