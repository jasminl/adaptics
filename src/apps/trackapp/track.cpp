#include "track.h"
#include "trackFeature.h"
#include "trackMatch.h"

#include <numeric>
#include <algorithm>
#include <iostream>

#ifdef _MANAGED
#pragma managed(push, off)
#endif

BOOL APIENTRY DllMain( HMODULE hModule,
                       DWORD  ul_reason_for_call,
                       LPVOID lpReserved
					 )
{
	switch (ul_reason_for_call)
	{
	case DLL_PROCESS_ATTACH:
	case DLL_THREAD_ATTACH:
	case DLL_THREAD_DETACH:
	case DLL_PROCESS_DETACH:
		break;
	}
    return TRUE;
}

#ifdef _MANAGED
#pragma managed(pop)
#endif

/***** trackMeanShift *****/

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
		m_model->expandScale(m_b, m_n);
		
		m_model->candidate(m_candidate,image,m_currentX,m_currentY,m_currentExpandedScale);		//Instantiation of candidate representation
		
		double rho0 = m_candidate->bhattacharyyaDistance(m_model,m_currentScale);	//Get initial rho
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

			pair<double,double> shift = m_candidate->spatialMeanShift();
			y1[0] += shift.first;
			y1[1] += shift.second;
			
			m_model->updateCandidate(m_candidate,image,m_currentX,m_currentY,m_currentExpandedScale);
			rho1 = m_candidate->bhattacharyyaDistance(m_model,m_currentScale);

			unsigned int meanShiftIter=1;
			while(rho1<rho0)
			{
				y1[0] = 0.5*(m_currentX + y1[0]);
				y1[1] = 0.5*(m_currentY + y1[1]);

				//Calculate p1 and correlate
				m_model->updateCandidate(m_candidate,image,m_currentX,m_currentY,m_currentExpandedScale);

				rho1 = m_candidate->bhattacharyyaDistance(m_model,m_currentScale);

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
			m_model->updateCandidate(m_candidate,image,m_currentX,m_currentY,m_currentExpandedScale);

			if(firstIteration==true)
			{
				//If it's the first iteration, set to initial position
				firstIteration = false;
				s1 = m_currentScale;
			}
			
			m_candidate->weight(m_model);
			double ds = m_candidate->scaleMeanShift();

			s1 = s1 * pow(m_b,ds);
			m_candidate->cScale() = s1;										//Rest scale to new location
			m_currentExpandedScale = m_candidate->expandScale(m_b,m_n);		//Adjust searched range

			m_model->updateCandidate(m_candidate,image,m_currentX,m_currentY,m_currentExpandedScale);
			double rho2 = m_candidate->bhattacharyyaDistance(m_model,m_candidate->currentScale());

			unsigned int meanShiftIter = 1;

			while(rho2<rho1)
			{
				s1 = 0.5 * (m_currentScale + s1);
				m_candidate->cScale() = s1;									//Rest scale to new location
				m_currentExpandedScale = m_candidate->expandScale(m_b,m_n); //Adjust searched range

				//Calculate p1 and correlate
				m_model->updateCandidate(m_candidate,image,m_currentX,m_currentY,m_currentExpandedScale);
				rho2 = m_candidate->bhattacharyyaDistance(m_model,m_candidate->currentScale());

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
			m_currentExpandedScale = m_candidate->expandScale(m_b,m_n);
		}

		
		nbIterAll++;
	}

	m_model->killCandidate(m_candidate);

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

/***** trackMSTargetBase  *****/
void trackMSTargetColor::weight(trackMSTargetBase* model)
{
	trackMSTargetColor* pModel = static_cast<trackMSTargetColor*>(model);

	pair<double*,unsigned int> p,s,r;
	pair<unsigned int*, unsigned int> q;

	for(unsigned int i=0;i<m_w2.size(); i++)
	{

		p = m_w2(i);
		q = m_bin2(i);
		s = m_hist(i); 
		r = pModel->m_hist(i);

		for(unsigned int j=0;j<p.second;j++)
		{
			if(s.first[q.first[j]] != 0)
			{
				p.first[j] = sqrt(r.first[q.first[j]]/s.first[q.first[j]]);
			}
			else
			{
				p.first[j] = 0;
			}
		}
	}
}

pair<double, double> trackMSTargetBase::spatialMeanShift()
{
	int n;
	vector<double> ss = setupScale(n);

	double numeratorX=0,numeratorY=0,denominator=0, xc, yc, skxwa, hsqn;
	vector<double>::iterator p,q,r,t;
	vector< pair<double,double> >::iterator u;

	for(p = m_s.begin(), q = ss.begin(); p != m_s.end(); p++, q++)
	{
		hsqn = hs(*q,n);
		if(hsqn == 0)
		{
			//The contribution of this scale (one of the two extremes) will be zero, so skip it
			continue;
		}
 
		pair<double*,unsigned int> r,s;
		pair<unsigned int*, unsigned int> q;
		r = m_dist2[*p];
		s = m_w2[*p];
		pair< pair<double,double>*,unsigned int> v = m_cenPix2[*p];

		for(unsigned int i=0;i<r.second;i++)
		{	
			skxwa = kp(r.first[i],*p) * s.first[i];	

			xc = skxwa * v.first[i].first;			//Here use unnormalized pixels coordinates (x)
			yc = skxwa * v.first[i].second;			//Here use unnormalized pixels coordinates (y)

			numeratorX	+= hsqn * xc;
			numeratorY	+= hsqn * yc;
			denominator += hsqn * abs(skxwa);
		}
	}

	return pair<double, double>(numeratorX/denominator,numeratorY/denominator);
}

double trackMSTargetBase::scaleMeanShift()
{
	int n;
	vector<double> ss = setupScale(n);

	double numerator =0, denominator = 0;

	vector<double>::iterator p,q,r,t;
	for(p = m_s.begin(), q = ss.begin(); p != m_s.end(); p++, q++)
	{
		double shxwa = 0;

		pair<double*,unsigned int> r = m_dist2[*p], s = m_w2[*p];

		for(unsigned int i=0;i<r.second;i++)
		{
			shxwa += hxq(r.first[i],*p) * s.first[i];
		}

		numerator += shxwa * *q;
		denominator += shxwa;
	}


	return numerator / denominator;
}

/***** trackMSTargetColor *****/

void trackMSTargetColor::allocateHistograms(vector<double> scale, bool backgroundWeight)
{

	//Calculate total number of bins (product of each dimension)
	unsigned int histogramSize = accumulate(m_nbBinsPerDim.begin(),m_nbBinsPerDim.end(),1,multiplies<unsigned int>());
	
	//For each scale, allocate one histogram
	unsigned int i;vector<double>::iterator p;
	for(p = scale.begin(), i=0; p != scale.end() ; p++, i++)
	{
		m_hist.changeIndex(i,*p);	//Change scale index
		m_hist.realloc(*p,histogramSize);
	}

	if(backgroundWeight == true)
	{
		//Allocate background histogram
		m_backHistogram = vector<double>(histogramSize);
	}
}


/***** trackMSTargetRGB *****/

trackMSTargetRGB::trackMSTargetRGB(double x, 
					  double y,
					  double hx,
					  double hy,
					  vector<double> s,
					  unsigned char* image,
					  unsigned int imageWidth,
					  unsigned int imageHeight,
					  double backSize,
					  unsigned int rNbBin,
					  unsigned int gNbBin,
					  unsigned int bNbBin):trackMSTargetColor(x,y,hx,hy,s,image,imageWidth,imageHeight,backSize)
{

	m_nbBinsPerDim.push_back(rNbBin);										//Setup histogram parameters
	m_nbBinsPerDim.push_back(gNbBin);
	m_nbBinsPerDim.push_back(bNbBin);
	
	m_nbLevelsPerBin.push_back(256.0 / m_nbBinsPerDim[0]);				//Get number of color levels per bin for each dimension
	m_nbLevelsPerBin.push_back(256.0 / m_nbBinsPerDim[1]);
	m_nbLevelsPerBin.push_back(256.0 / m_nbBinsPerDim[2]);

	allocateHistograms(s,(backSize>0)?true:false);					//Allocate histogram space

	makeHistogram(m_s);												//Make histograms
}

double trackMSTargetColor::bhattacharyyaDistance(trackMSTargetBase* a, double scale)
{
	trackMSTargetColor* pa = dynamic_cast<trackMSTargetColor*>(a);

	//Obtain histogram of 'a', using the fact it is always a model, and that the scale 1 is important
	pair<double*,unsigned int> r = pa->m_hist.centerScale();
	
	//map<double, vector<double> >::iterator t = m_histogram.find(scale);			//Obtain 'current' histogram
	pair<double*,unsigned int> t = m_hist[scale];
		

	double sum = 0;

	for(unsigned int i=0;i<t.second;i++)
	{
		sum += sqrt(t.first[i] * r.first[i]);
	}
	
	return (sum>1?1:sum);
}

void trackMSTargetColor::makeHistogram(vector<double>& scale)
{
	for(vector<double>::iterator ps = scale.begin(); ps != scale.end() ; ps++)
	{
		//Find histogram that corresponds to current scale
		double* q = m_hist[*ps].first;

		//Calculate target ellipse size
		double majAxis = m_hy * *ps / 2;							//Vertical axis
		double minAxis = m_hx * *ps / 2;							//Horizontal axis
		double majAxis2 = majAxis * majAxis;
		double minAxis2 = minAxis * minAxis;

		//Background ellipse
		double backMajAxis  = m_backSize * majAxis;
		double backMinAxis  = m_backSize * minAxis;
		double backMajAxis2 = backMajAxis * backMajAxis;
		double backMinAxis2 = backMinAxis * backMinAxis;
		
		unsigned int nbBytesPerPosition = (unsigned)m_nbBinsPerDim.size();	//This tells us how many dimensions the image has

		double backSum=0;													//Total nb of pixels in background histogram
		
		if (m_backSize > 0)
		{
			//Use background weighting
			boundingRect(backMajAxis,backMinAxis);							//Coordinates of rectangle must include background area
 			m_backHistogram.assign(m_backHistogram.size(),0);				//Reset background histogram to zeros
		}
		else
		{
			//No background weighting
			boundingRect(majAxis,minAxis);
		}
 		
		unsigned char* pImage = m_image + nbBytesPerPosition * (m_vertical1 * m_width + m_horizontal1);	//Advance position to beginning of bounding rectangle

		unsigned int* pb = m_bin2.realloc(*ps,static_cast<unsigned int>(ceil(2*majAxis * 2*minAxis)));
		pair<double,double>* pcp= m_cenPix2.realloc(*ps,static_cast<unsigned int>(ceil(2*majAxis * 2*minAxis)));
		double* pd = m_dist2.realloc(*ps,static_cast<unsigned int>(ceil(2*majAxis * 2*minAxis)));

		unsigned int i,j,bin;
		double nx,ny,dist,kx;

		int nb=0;
		unsigned int newSize=0;
		//First check it's in the background, and then whether it's in the target area
		for(i=m_vertical1;i<=m_vertical2;i++)
		{
			for(j=m_horizontal1;j<=m_horizontal2;j++, pImage+=nbBytesPerPosition)
			{
				//Check that pixel belongs to background and or ellipse
				if ( (i-m_y)*(i-m_y)*majAxis2 + (j-m_x)*(j-m_x)*minAxis2  <= majAxis2*minAxis2)
				{
					//Add to target histogram
						
					pcp->first  = j-m_x;									//Stored centered, unnormalized pixels
					pcp->second = i-m_y;
					pcp++;
					
					nx = (j-m_x) / m_hx;									//Normalized and centered coordinates
					ny = (i-m_y) / m_hy;
					
					dist = sqrt(nx*nx + ny*ny);								//Euclidean distance

					*pd++ = dist;											//Save distance

					kx = kp(dist,*ps);										//Kernel value at current scale
					
					bin = quantize(*pImage,*(pImage+1),*(pImage+2));		//Histogram bin

					*pb++ = bin;											//Save histogram bin of current pixel
					
					q[bin] += kx;

					newSize++;
				}
				else if( (m_backSize>0) && ( (i-m_y)*(i-m_y)*backMajAxis2 + (j-m_x)*(j-m_x)*backMinAxis2  <= backMajAxis2*backMinAxis2) )
				{
					//Add to background histogram
					bin = quantize(*pImage,*(pImage+1),*(pImage+2));	//Histogram bin

					m_backHistogram[bin] += 1;							//Add to target histogram
					backSum += 1;										//Add to normalization factor
				}
			}

			//Shift by appropriate horizontal offset
			pImage += nbBytesPerPosition * ( m_width - (m_horizontal2 - m_horizontal1) - 1);
		}
	
		//Resize m_bin given the actual number of pixels in the target ellipse, and allocate m_w accordingly
		m_bin2.actualSize(*ps) = newSize;
		m_w2.realloc(*ps,newSize);
		m_cenPix2.actualSize(*ps) = newSize;
		m_dist2.actualSize(*ps) = newSize;
	 
		if(m_backSize>0 && backSum!=0)
		{
			double minValue = 10000000;

			//Adjust background histogram
			vector<double>::iterator p;
			for(p = m_backHistogram.begin() ; p != m_backHistogram.end() ; p++)
			{
				*p /= backSum;											//Normalize such that sum is 1

				if(*p > 0 && *p<minValue)								//Find smallest nonzero element
				{
					minValue = *p;
				}
			}

			vector<double> vu(m_backHistogram.size(),1);						//'vu' quantity for background weighting
			vector<double>::iterator r = vu.begin();

			for(p=m_backHistogram.begin();p!=m_backHistogram.end();p++,r++)
			{
				//Get weights similar to backprojection
				if(*p > 0 )
				{
					*r = min(minValue / *p, 1);							//Note, if *p=0, automatically vu=1
				}
			}

			//Background weighting
			pImage = m_image + nbBytesPerPosition * (m_vertical1 * m_width + m_horizontal1);		//Reset position
			
			double c = 0;																			//Normalization constant
			
			boundingRect(majAxis,minAxis);															//Always use the target area as bounding rectangle

			for(i=m_vertical1;i<=m_vertical2;i++)
			{
				for(j=m_horizontal1;j<=m_horizontal2;j++, pImage+=nbBytesPerPosition)
				{
					//Check that pixel belongs to background and or ellipse
					if ( (i-m_y)*(i-m_y)*majAxis2 + (j-m_x)*(j-m_x)*minAxis2  <= majAxis2*minAxis2)
					{
						//Add to target histogram
						
						double nx = (j-m_x) / m_hx;							//Normalized and centered coordinates
						double ny = (i-m_y) / m_hy;
						
						double dist = sqrt(nx*nx + ny*ny);					//Euclidean distance

						double kx = kp(dist,*ps);							//Kernel value at current scale
						
						bin = quantize(*pImage,*(pImage+1),*(pImage+2));	//Histogram bin

						c += kx * vu[bin];									//Add to normalization constant
					}
				}

				//Shift by appropriate horizontal offset
				pImage += nbBytesPerPosition * ( m_width - (m_horizontal2 - m_horizontal1) - 1);
			}

			//Adjusted previously calculated target histogram
			unsigned int i;
			for(i=0, r=vu.begin() ; r<vu.end() ; i++,r++)
			{
				q[i] *= *r/c;
			}
		}
	}
}

/***** Misc functions *****/
vector<bool> track1Frame(pair<double,double>& location,
				 double& scale,
				 unsigned char* image,
				 vector<trackFlow*> tracker,
				 double hx,
   			     double hy,
				 unsigned int width,
				 unsigned int height)
{
	//Create boolean vector indicating whether the tracked object is valid, for each tracking object
 	vector<bool> validTrack(tracker.size());
	vector<bool>::iterator q = validTrack.begin();

	double denominator=0,numerator=0;				//Qties used to average output 
	location.first = location.second = scale = 0;	//Initialize output variables to zero

	for(vector<trackFlow*>::iterator p=tracker.begin() ; p != tracker.end() ; p++, q++)	
	{
		*q = (*p)->track(image);

		if(*q == true)
		{
			//If this a valid target, include it in the fused x,y,scale estimates.
			trackFlow::coord current = (*p)->currentCoordinates();		
			location.first += current.s_x;
			location.second += current.s_y;
			scale += current.s_scale;

			denominator++;
		}
	}

	//Average output
	location.first /= denominator;
	location.second /= denominator;
	scale /= denominator;

	//Create cropped image for shape matching
	unsigned char* cropped = NULL;
	
	crop(cropped,scale*hx,scale*hy,location.first,location.second,image,width,height);
	
	//Perform shape matching if applicable
	for(vector<trackFlow*>::iterator p=tracker.begin() ; p != tracker.end() ; p++)	
	{
		if((*p)->featFilter() == NULL) continue;		//Skip this step if no feature descriptor has been defined

		//Use target filter to compute and process keypoints
		(*p)->targetFilter()->compute(cropped,scale*hx,scale*hy);
	//	(*p)->targetFilter()->process();
		
		//DEBUGGIN
		(*p)->featFilter()->feature()->show();		//THIS SHOWS that zero points are in the array
		(*p)->targetFilter()->feature()->show();

TODO: show features on image
		//Perform matching between target and model
		(*p)->matcher()->match((*p)->featFilter()->feature(),(*p)->targetFilter()->feature());

TODO: show matches on image
	}

	delete[] cropped;	//Clean up matching image

	return validTrack;	//Return vector indicating validity of tracked object
}

unsigned char* crop(unsigned char*& buffer, unsigned int width, unsigned int height,double x,double y,unsigned char* image,unsigned int imageWidth,unsigned int imageHeight)
{
	if (buffer != NULL)
	{
		delete[] buffer;
	}

	//Determine coordinates of area to crop
	unsigned int x1 = max(0,x-width/2);	
	unsigned int x2 = min(imageWidth,x+width/2);
	unsigned int y1 = max(0,y-height/2);
	unsigned int y2 = min(imageHeight,y+height/2);

	//Calculate range
	double xRange = x2 - x1 + 1;
	double yRange = y2 - y1 + 1;

	//Allocate buffer
	buffer = new unsigned char[3 * (unsigned int)(ceil(xRange) * ceil(yRange))];

	//Copy image
	unsigned char* p = buffer;								//Point to cropped image
	unsigned char* q = image + y1 * imageWidth + x1;		//Pointer to base image
	for(unsigned int i=y1;i<y2;i++)
	{
		memcpy(p,q,3 * xRange * sizeof(unsigned char));		//Copy 1 row
		q += 3 * imageWidth;								//Forward to next row
	}

	return buffer;											//Return result
}
