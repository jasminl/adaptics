#include "trackMSTargetColor.h"

#include <cmath>
#include <algorithm>

using namespace std;

pair<vector<unsigned int>, vector<double> > trackMSTargetColor::histParams() const
{
	return pair<vector<unsigned int> , vector<double> >(m_nbBinsPerDim,m_nbLevelsPerBin);
}

void trackMSTargetColor::setHistogram(vector<unsigned int> nb, vector<double> sz)
{
	m_nbBinsPerDim = nb;
	m_nbLevelsPerBin = sz;
}

double trackMSTargetColor::bhattacharyya_distance(TrackMSTargetBase* a, double scale)
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
		double majAxis = _hy * *ps / 2;							//Vertical axis
		double minAxis = _hx * *ps / 2;							//Horizontal axis
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
			bounding_rect(backMajAxis,backMinAxis);							//Coordinates of rectangle must include background area
 			m_backHistogram.assign(m_backHistogram.size(),0);				//Reset background histogram to zeros
		}
		else
		{
			//No background weighting
			bounding_rect(majAxis,minAxis);
		}

		unsigned char* pImage = _image + nbBytesPerPosition * (_vertical1 * _width + _horizontal1);	//Advance position to beginning of bounding rectangle

		unsigned int* pb = _bin2.realloc(*ps,static_cast<unsigned int>(ceil(2*majAxis * 2*minAxis)));
		pair<double,double>* pcp= _cenPix2.realloc(*ps,static_cast<unsigned int>(ceil(2*majAxis * 2*minAxis)));
		double* pd = _dist2.realloc(*ps,static_cast<unsigned int>(ceil(2*majAxis * 2*minAxis)));

		unsigned int i,j,bin;
		double nx,ny,dist,kx;

		int nb=0;
		unsigned int newSize=0;
		//First check it's in the background, and then whether it's in the target area
		for(i=_vertical1;i<=_vertical2;i++)
		{
			for(j=_horizontal1;j<=_horizontal2;j++, pImage+=nbBytesPerPosition)
			{
				//Check that pixel belongs to background and or ellipse
				if ( (i-_y)*(i-_y)*majAxis2 + (j-_x)*(j-_x)*minAxis2  <= majAxis2*minAxis2)
				{
					//Add to target histogram

					pcp->first  = j-_x;									//Stored centered, unnormalized pixels
					pcp->second = i-_y;
					pcp++;

					nx = (j-_x) / _hx;									//Normalized and centered coordinates
					ny = (i-_y) / _hy;

					dist = sqrt(nx*nx + ny*ny);								//Euclidean distance

					*pd++ = dist;											//Save distance

					kx = kp(dist,*ps);										//Kernel value at current scale

					bin = quantize(*pImage,*(pImage+1),*(pImage+2));		//Histogram bin

					*pb++ = bin;											//Save histogram bin of current pixel

					q[bin] += kx;

					newSize++;
				}
				else if( (m_backSize>0) && ( (i-_y)*(i-_y)*backMajAxis2 + (j-_x)*(j-_x)*backMinAxis2  <= backMajAxis2*backMinAxis2) )
				{
					//Add to background histogram
					bin = quantize(*pImage,*(pImage+1),*(pImage+2));	//Histogram bin

					m_backHistogram[bin] += 1;							//Add to target histogram
					backSum += 1;										//Add to normalization factor
				}
			}

			//Shift by appropriate horizontal offset
			pImage += nbBytesPerPosition * ( _width - (_horizontal2 - _horizontal1) - 1);
		}

		//Resize m_bin given the actual number of pixels in the target ellipse, and allocate m_w accordingly
		_bin2.actualSize(*ps) = newSize;
		_w2.realloc(*ps,newSize);
		_cenPix2.actualSize(*ps) = newSize;
		_dist2.actualSize(*ps) = newSize;

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
					*r = min(minValue / *p, 1.0);							//Note, if *p=0, automatically vu=1
				}
			}

			//Background weighting
			pImage = _image + nbBytesPerPosition * (_vertical1 * _width + _horizontal1);		//Reset position

			double c = 0;																			//Normalization constant

			bounding_rect(majAxis,minAxis);															//Always use the target area as bounding rectangle

			for(i=_vertical1;i<=_vertical2;i++)
			{
				for(j=_horizontal1;j<=_horizontal2;j++, pImage+=nbBytesPerPosition)
				{
					//Check that pixel belongs to background and or ellipse
					if ( (i-_y)*(i-_y)*majAxis2 + (j-_x)*(j-_x)*minAxis2  <= majAxis2*minAxis2)
					{
						//Add to target histogram

						double nx = (j-_x) / _hx;							//Normalized and centered coordinates
						double ny = (i-_y) / _hy;

						double dist = sqrt(nx*nx + ny*ny);					//Euclidean distance

						double kx = kp(dist,*ps);							//Kernel value at current scale

						bin = quantize(*pImage,*(pImage+1),*(pImage+2));	//Histogram bin

						c += kx * vu[bin];									//Add to normalization constant
					}
				}

				//Shift by appropriate horizontal offset
				pImage += nbBytesPerPosition * ( _width - (_horizontal2 - _horizontal1) - 1);
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

void trackMSTargetColor::weight(TrackMSTargetBase* model)
{
	trackMSTargetColor* pModel = static_cast<trackMSTargetColor*>(model);

	pair<double*,unsigned int> p,s,r;
	pair<unsigned int*, unsigned int> q;

	for(unsigned int i=0;i<_w2.size(); i++)
	{

		p = _w2(i);
		q = _bin2(i);
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


trackMSTargetColor::trackMSTargetColor(double x,
				  double y,
				  double hx,
				  double hy,
				  std::vector<double> s,
				  unsigned char* image,
				  unsigned int imageWidth,
				  unsigned int imageHeight,
				  double backSize):TrackMSTargetBase(x,y,hx,hy,s,image,imageWidth,imageHeight),m_backSize(backSize)
{}
