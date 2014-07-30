#include "TrackMSTargetBase.h"
#include <cmath>

using namespace std;

#define PI 3.1416

TrackMSTargetBase::TrackMSTargetBase()
: _x(0), _y(0), _hx(0), _hy(0),
  _s(vector<double>()), _image(nullptr), _width(0), _height(0),
  _vertical1(0), _vertical2(0), _horizontal1(0), _horizontal2(0)
{};

TrackMSTargetBase::TrackMSTargetBase(double x, double y, double hx, double hy,
				  std::vector<double> s, unsigned char* image, unsigned int imageWidth, unsigned int imageHeight)
: _x(x), _y(y), _hx(hx), _hy(hy),
  _s(s), _image(image), _width(imageWidth), _height(imageHeight),
  _vertical1(0), _vertical2(0), _horizontal1(0), _horizontal2(0)
{
	allocate(s);
};

void TrackMSTargetBase::bounding_rect(double verticalAxis, double horizontalAxis)
{
	_vertical1   = static_cast<unsigned int>(max(_y - verticalAxis,0.0));
	_vertical2   = static_cast<unsigned int>(min(_y + verticalAxis,_height-1.0));
	_horizontal1 = static_cast<unsigned int>(max(_x - horizontalAxis,0.0));
	_horizontal2 = static_cast<unsigned int>(min(_x + horizontalAxis,_width-1.0));
}

pair<double, double> TrackMSTargetBase::spatial_meanshift()
{
	int n;
	vector<double> ss = setup_scale(n);

	double numeratorX = 0, numeratorY = 0, denominator = 0;
	for(auto p = _s.begin(), q = ss.begin(); p != _s.end(); p++, q++)
	{
		double hsqn = hs(*q, n);
		if(hsqn == 0) //The contribution of this scale (one of the two extremes) will be zero, so skip it
			continue;

		auto r = _dist2[*p];
		auto s = _w2[*p];
		auto v = _cenPix2[*p];

		for(unsigned int i = 0;i < r.second; i++)
		{
			double skxwa = kp(r.first[i], *p) * s.first[i];

			double xc = skxwa * v.first[i].first;	//Here use unnormalized pixels coordinates (x)
			double yc = skxwa * v.first[i].second;	//Here use unnormalized pixels coordinates (y)

			numeratorX	+= hsqn * xc;
			numeratorY	+= hsqn * yc;
			denominator += hsqn * abs(skxwa);
		}
	}
	return {numeratorX/denominator, numeratorY/denominator};
}

double TrackMSTargetBase::scale_meanshift()
{
	int n;
	vector<double> ss = setup_scale(n);

	double numerator = 0, denominator = 0;
	for(auto p = _s.begin(), q = ss.begin(); p != _s.end(); p++, q++)
	{
		double shxwa = 0;

		auto r = _dist2[*p], s = _w2[*p];
		for(unsigned int i = 0; i < r.second; i++)
			shxwa += hxq(r.first[i],*p) * s.first[i];
		numerator += shxwa * *q;
		denominator += shxwa;
	}
	return numerator / denominator;
}

void TrackMSTargetBase::allocate(std::vector<double>& s)
{
	for(unsigned int i = 0;i < _bin2.size(); i++)
	{
		_bin2.changeIndex(i, s[i]);
		_w2.changeIndex(i, s[i]);
		_cenPix2.changeIndex(i, s[i]);
		_dist2.changeIndex(i, s[i]);
	}
}

double TrackMSTargetBase::current_scale() const
{
	unsigned int index = static_cast<unsigned int>(_s.size()/2);
	return _s[index];
}

void TrackMSTargetBase::set(double x, double y, const std::vector<double> s)
{
	_x = x;
	_y = y;
	_s = s;
}

double& TrackMSTargetBase::c_scale()
{
	unsigned int index = static_cast<unsigned int>(_s.size()/2);
	return _s[index];
}

vector<double> TrackMSTargetBase::expand_scale(double b, int n)
{
	double c = current_scale();		//Obtain current scale (if there is only one entry, returns it)

	_s = vector<double>(2*n+1);	//Allocate space for new scale std::vector
	n = -n;
	for(auto p = _s.begin(); p != _s.end(); p++, n++)
		*p = c * pow(b, n);
	return _s;
}

double TrackMSTargetBase::epanichnikov(double distance, double cd, double d)
{
	if(distance <= 1)
		return 0.5 / cd * (d + 2)*(1 - distance);
	return 0;
}

double TrackMSTargetBase::hs(double scale, double n)
{
	if (n == 0)
		return 1;
	return 1 - (scale/n)*(scale/n);
}

double TrackMSTargetBase::kp(double distance, double sigma)
{
	double s1 = sigma / 1.6;
	double s2 = sigma * 1.6;

	double term1 = exp((- distance * distance)/(2 * s1 * s1)) / (2 * PI * pow(s1,4.0));
	double term2 = exp((- distance * distance)/(2 * s2 * s2)) / (2 * PI * pow(s2,4.0));

	return term1 - term2;
}

double TrackMSTargetBase::hxq(double distance, double sigma)
{
	double fd = 2*PI*sigma*sigma / 1.6;	//Scale divided by 1.6
	double ft = 2*PI*sigma*sigma * 1.6;	//Scale times 1.6

	double fdex = 2*sigma*sigma / 1.6;	//Scale divided by 1.6
	double ftex = 2*sigma*sigma * 1.6;	//Scale times 1.6

	double term1 = exp(-distance*distance / fdex ) / fd;
	double term2 = exp(-distance*distance / ftex ) / ft;

	return term1 - term2;
}

double TrackMSTargetBase::bhattacharyya_distance(TrackMSTargetBase* a, double scale)
{
	return 0;
}

vector<double> TrackMSTargetBase::setup_scale(int& n)
{
	vector<double> ss;

	if(_s.size() == 1)
	{
		//Mean shift is only done in space, not scale
		ss.push_back(0);
		n = 0;
	}
	else
	{
		//Mean shift is performed in scale space
		n  = static_cast<int>(_s.size() / 2);	//Get radius of tracked scales
		ss = vector<double>(_s.size());

		double pn = -n;
		for(auto p = ss.begin(); p != ss.end(); p++)
		{
			*p = pn++;
		}
	}

	return ss;
}

vector<double>& TrackMSTargetBase::scale_range()
{
	return _s;
}

pair<unsigned int, unsigned int> TrackMSTargetBase::wh() const
{
	return std::pair<unsigned int,unsigned int>(_width,_height);
}

