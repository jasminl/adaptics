#include "TrackMSTargetColor.h"

#include <cmath>
#include <algorithm>

using namespace std;

TrackMSTargetColor::TrackMSTargetColor(double x, double y, double hx, double hy,
		vector<double> s, unsigned char* image, int image_width, int image_height,
		double back_size)
: TrackMSTargetBase(x, y, hx, hy, s, image, image_width, image_height),
  _back_size(back_size)
{}

pair<vector<unsigned int>, vector<double>> TrackMSTargetColor::histParams() const
{
	return make_pair(_bins_per_dim, _levels_per_bin);
}

void TrackMSTargetColor::set_histogram(vector<unsigned int> nb, vector<double> sz)
{
	_bins_per_dim = nb;
	_levels_per_bin = sz;
}

double TrackMSTargetColor::bhattacharyya_distance(TrackMSTargetBase* a, double scale)
{
	TrackMSTargetColor* pa = dynamic_cast<TrackMSTargetColor*>(a);

	/*
	 * Obtain histogram of 'a', using the fact it is always a model,
	 * and that the scale 1 is important
	 */
	pair<double*,unsigned int> r = pa->_hist.centerScale();

	pair<double*,unsigned int> t = _hist[scale];

	double sum = 0;
	for(unsigned int i = 0; i < t.second; i++)
		sum += sqrt(t.first[i] * r.first[i]);

	return (sum > 1?1:sum);
}

void TrackMSTargetColor::makeHistogram(vector<double>& scale)
{
	for(auto ps = scale.begin(); ps != scale.end() ; ps++)
	{
		//Find histogram that corresponds to current scale
		double* q = _hist[*ps].first;

		//Calculate target ellipse size
		double majAxis = _hy * *ps / 2;							//Vertical axis
		double minAxis = _hx * *ps / 2;							//Horizontal axis
		double majAxis2 = majAxis * majAxis;
		double minAxis2 = minAxis * minAxis;

		//Background ellipse
		double backMajAxis  = _back_size * majAxis;
		double backMinAxis  = _back_size * minAxis;
		double backMajAxis2 = backMajAxis * backMajAxis;
		double backMinAxis2 = backMinAxis * backMinAxis;

		unsigned int nbBytesPerPosition = (unsigned)_bins_per_dim.size();	//This tells us how many dimensions the image has

		double backSum=0;													//Total nb of pixels in background histogram

		if (_back_size > 0)
		{
			//Use background weighting
			bounding_rect(backMajAxis, backMinAxis); //Coordinates of rectangle must include background area
 			_back_hist.assign(_back_hist.size(), 0); //Reset background histogram to zeros
		}
		else //No background weighting
			bounding_rect(majAxis,minAxis);

		unsigned char* pImage = _image + nbBytesPerPosition * (_vertical1 * _width + _horizontal1);	//Advance position to beginning of bounding rectangle

		unsigned int* pb = _bin2.realloc(*ps, static_cast<unsigned int>(ceil(2 * majAxis * 2 * minAxis)));
		auto pcp = _cenPix2.realloc(*ps, static_cast<unsigned int>(ceil(2 * majAxis * 2 * minAxis)));
		double* pd = _dist2.realloc(*ps, static_cast<unsigned int>(ceil(2 * majAxis * 2 * minAxis)));

		unsigned int i,j,bin;
		double nx,ny,dist,kx;

		unsigned int new_size = 0;
		//First check it's in the background, and then whether it's in the target area
		for(i = _vertical1; i <= _vertical2; i++)
		{
			for(j = _horizontal1; j <= _horizontal2; j++, pImage += nbBytesPerPosition)
			{
				//Check that pixel belongs to background and or ellipse
				if ( (i-_y)*(i-_y)*majAxis2 + (j-_x)*(j-_x)*minAxis2  <= majAxis2*minAxis2)
				{
					//Add to target histogram
					pcp->first  = j - _x;		//Stored centered, unnormalized pixels
					pcp->second = i - _y;
					pcp++;

					nx = (j-_x) / _hx;  		//Normalized and centered coordinates
					ny = (i-_y) / _hy;

					dist = sqrt(nx*nx + ny*ny);		//Euclidean distance

					*pd++ = dist;					//Save distance

					kx = kp(dist,*ps); 	//Kernel value at current scale

					bin = quantize(*pImage, *(pImage+1), *(pImage+2));		//Histogram bin

					*pb++ = bin;	//Save histogram bin of current pixel

					q[bin] += kx;

					new_size++;
				}
				else if( (_back_size>0) && ( (i-_y)*(i-_y)*backMajAxis2 + (j-_x)*(j-_x)*backMinAxis2  <= backMajAxis2*backMinAxis2) )
				{
					//Add to background histogram
					bin = quantize(*pImage,*(pImage+1),*(pImage+2));	//Histogram bin

					_back_hist[bin] += 1;	//Add to target histogram
					backSum += 1;	//Add to normalization factor
				}
			}

			//Shift by appropriate horizontal offset
			pImage += nbBytesPerPosition * ( _width - (_horizontal2 - _horizontal1) - 1);
		}

		//Resize m_bin given the actual number of pixels in the target ellipse, and allocate m_w accordingly
		_bin2.actual_size(*ps) = new_size;
		_w2.realloc(*ps,new_size);
		_cenPix2.actual_size(*ps) = new_size;
		_dist2.actual_size(*ps) = new_size;

		if(_back_size>0 && backSum!=0)
		{
			double minValue = 10000000;

			//Adjust background histogram
			//vector<double>::iterator p;
			for(auto p = _back_hist.begin() ; p != _back_hist.end() ; p++)
			{
				*p /= backSum;											//Normalize such that sum is 1

				if(*p > 0 && *p<minValue)								//Find smallest nonzero element
				{
					minValue = *p;
				}
			}

			vector<double> vu(_back_hist.size(),1);						//'vu' quantity for background weighting
			vector<double>::iterator r = vu.begin();

			for(auto p=_back_hist.begin();p!=_back_hist.end();p++,r++)
			{
				//Get weights similar to backprojection
				if(*p > 0 )
					*r = min(minValue / *p, 1.0);							//Note, if *p=0, automatically vu=1
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

void TrackMSTargetColor::weight(TrackMSTargetBase* model)
{
	TrackMSTargetColor* pModel = static_cast<TrackMSTargetColor*>(model);
	for(unsigned int i = 0; i < _w2.size(); i++)
	{

		auto p = _w2(i);
		auto q = _bin2(i);
		auto s = _hist(i);
		auto r = pModel->_hist(i);

		for(unsigned int j = 0; j <p.second; j++)
		{
			if(s.first[q.first[j]] != 0)
				p.first[j] = sqrt(r.first[q.first[j]]/s.first[q.first[j]]);
			else
				p.first[j] = 0;
		}
	}
}

void TrackMSTargetColor::allocateHistograms(vector<double> scale, bool backgroundWeight)
{

	//Calculate total number of bins (product of each dimension)
	unsigned int histogramSize = accumulate(_bins_per_dim.begin(),_bins_per_dim.end(),
			1, multiplies<unsigned int>());

	//For each scale, allocate one histogram
	unsigned int i = 0;
	for(auto p: scale)
	{	//Change scale index
		_hist.changeIndex(i,p);
		_hist.realloc(p, histogramSize);
		 i++;
	}


	if(backgroundWeight == true)	//Allocate background histogram
		_back_hist = vector<double>(histogramSize);
}
