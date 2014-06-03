#include "visio/Layer.h"
#include "tinyxml/tinyxml.h"
#include <omp.h>
#include <time.h>
#include <fstream>
#include <stdexcept>

using namespace std;

namespace visio
{
CLayer::CLayer(unsigned int VSize, unsigned int HSize, unsigned int totFrames, unsigned int nbDirections,
		string name, bool record, Buffer* buffer)
: _vSize(VSize), _hSize(HSize), _nb_frames(totFrames), _nb_dir(nbDirections),
  _level_name(name), _record(record), _buf(buffer)
{
	if(_nb_dir % 2)
		throw runtime_error("CLayer::CLayer: uneven number of directions");
}

unsigned int CLayer::chop_input2(unsigned int vsize, unsigned int hsize)
{
	unsigned int totSize = vsize*hsize;

	if(totSize%2 == 0)
	{
		//Even sized input
		return totSize/2;
	}
	else
	{
		//Odd sized input: in this case the first half is one unit less than the second half
		return (totSize-1)/2;
	}
}

unsigned int CLayer::chop_input2(unsigned int totSize)
{
	if(totSize%2 == 0)
	{
		//Even sized input
		return totSize/2;
	}
	else
	{
		//Odd sized input: in this case the first half is one unit less than the second half
		return (totSize-1)/2;
	}
}

void CLayer::show_activity(double* array) const
{
	for(unsigned int i=0;i<_vSize;i++)
	{
		for(unsigned int j=0;j<_hSize;j++)
		{
			cout<<array[i*_hSize+j]<<" ";
		}
		cout<<endl;
	}

}

double* CLayer::create_array(double*& pointer, unsigned int vSize, unsigned int hSize)
{
	pointer = new double[vSize*hSize];
	memset(pointer,0,sizeof(double)*vSize*hSize);
	return pointer;
}

void CLayer::delete_array(double*& pointer)
{
	delete[] pointer;
}

void CLayer::set_value(double* array, double value)
{
	for(unsigned int i=0;i<_vSize;i++)
	{
		for(unsigned int j=0;j<_hSize;j++)
		{
			array[i*_hSize+j] = value;
		}
	}
}
}
