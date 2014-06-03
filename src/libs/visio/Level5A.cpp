#include "visio/Level5A.h"
#include "visio/Level4.h"
#include <stdexcept>
#include <omp.h>

using namespace std;

namespace visio
{
CLevel5A::CLevel5A(int iter_per_frame, double dt,
			double A7, double Ke, double Kb, double Kz,
			unsigned int verticalSize, unsigned int horizontalSize,
			unsigned int totFrames,
			unsigned int nbDirections,
		    double sigma,
			double I,
			unsigned int nbV2Depth,
			Input<int>* inputLayer, const std::string& output_folder,
			string name, bool record, Buffer* buffer, bool testKernels)
: KernelLayer(verticalSize, horizontalSize, totFrames, nbDirections, sigma, sigma,
		sigma, sigma, name, output_folder, iter_per_frame, dt, record, buffer), _nb_depth(nbV2Depth)
{
	_cst21 = 1-dt*A7;
	_cst22 = dt*A7*Ke;
	_cst23 = dt*A7*Kz;
	_cst24 = -dt*A7*Kb;

	//Weight filter
	_d0s1K.weight(I / (2*PI*sigma*sigma) );
	_d0s1K.set_intrinsic_weight(I / (2*PI*sigma*sigma) );

	//Delete redundant filters (only keep m_d0s1K)
	_d1s1K.remove();
	_d2s1K.remove();
	_d3s1K.remove();

	_d0s2K.remove();
	_d1s2K.remove();
	_d2s2K.remove();
	_d3s2K.remove();

#if NUM_THREADS == 4 || NUM_THREADS == 16
	_d5s1K.remove();
	_d7s1K.remove();
	_d5s2K.remove();
	_d7s2K.remove();
#endif

	//SVD for horizontal and vertical kernels
	_d0s1K.svd();

	//Allocate V2 boundary buffer space
	create_array(_V2, nbV2Depth, _hSize * _vSize);

	//Assign input layer
	_input = inputLayer;

	setBaseB();
	setBaseSB();

	if(testKernels)
		time_kernel(verticalSize, horizontalSize);
}

CLevel5A::~CLevel5A()
{
	delete_array(_V2,_nb_depth);

	for(auto p = _input_file.begin() ; p != _input_file.end() ; p++)
	{
		if((*p)->is_open())
			(*p)->close();
	}
}


void CLevel5A::setBaseB()
{
#if NUM_THREADS == 1
	_baseB.push_back(_buf->_buffer0);
#elif NUM_THREADS == 4
	_baseB.push_back(_buf->_buffer0);
	_baseB.push_back(_buf->_buffer1);
	_baseB.push_back(_buf->_buffer2);
	_baseB.push_back(_buf->_buffer3);
#elif NUM_THREADS == 16
	_baseB.push_back(_buf->_buffer0);
	_baseB.push_back(_buf->_buffer1);
	_baseB.push_back(_buf->_buffer2);
	_baseB.push_back(_buf->_buffer3);
#endif
}
void CLevel5A::setBaseSB()
{
#if NUM_THREADS == 1
	_baseSB.push_back(_buf->_buffer1);
#elif NUM_THREADS == 4
	_baseSB.push_back(_buf->_buffer4);
	_baseSB.push_back(_buf->_buffer5);
#elif NUM_THREADS == 16
	_baseSB.push_back(_buf->_buffer4);
	_baseSB.push_back(_buf->_buffer5);
#endif
}

void CLevel5A::compute(void* object)
{
	CLevel4* L4 = static_cast<CLevel4*>(object);

	//Load V2 boundaries
	_input->V2(1,_V2[0]);
	_input->V2(2,_V2[1]);

#if NUM_THREADS == 1

	//Scale 1
	_d0s1K.conv_sep(_V2[0],_vSize,_hSize,_baseB[0],_baseSB[0]);
	update(m_d0s1,L4->array("D0S1"),_V2[0],_baseB[0]);
	update(m_d1s1,L4->array("D1S1"),_V2[0],_baseB[0]);
	update(m_d2s1,L4->array("D2S1"),_V2[0],_baseB[0]);
	update(m_d3s1,L4->array("D3S1"),_V2[0],_baseB[0]);
	update(m_d4s1,L4->array("D4S1"),_V2[0],_baseB[0]);
	update(m_d5s1,L4->array("D5S1"),_V2[0],_baseB[0]);
	update(m_d6s1,L4->array("D6S1"),_V2[0],_baseB[0]);
	update(m_d7s1,L4->array("D7S1"),_V2[0],_baseB[0]);

	//Scale 2
	_d0s1K.conv_sep(_V2[1],_vSize,_hSize,_baseB[0],_baseSB[0]);

	update(m_d0s2,L4->array("D0S2"),_V2[1],_baseB[0]);
	update(m_d1s2,L4->array("D1S2"),_V2[1],_baseB[0]);
	update(m_d2s2,L4->array("D2S2"),_V2[1],_baseB[0]);
	update(m_d3s2,L4->array("D3S2"),_V2[1],_baseB[0]);
	update(m_d4s2,L4->array("D4S2"),_V2[1],_baseB[0]);
	update(m_d5s2,L4->array("D5S2"),_V2[1],_baseB[0]);
	update(m_d6s2,L4->array("D6S2"),_V2[1],_baseB[0]);
	update(m_d7s2,L4->array("D7S2"),_V2[1],_baseB[0]);

#elif NUM_THREADS == 4 || NUM_THREADS == 16

	vector<double*> a = _a;
	vector<double*> pb = _baseB;
	vector<double*> psb = _baseSB;

	double* V2D0 = _V2[0];
	double* V2D1 = _V2[1];

	//First do convolutions
	omp_set_num_threads(2);

	Kernel<double>* kern = &_d0s1K;
	double* v2_0 = _V2[0];
	double* v2_1 = _V2[1];

	#pragma omp parallel shared(kern, pb, psb)
	{
		int id = omp_get_thread_num();
		switch(id)
		{
		case 0:
			kern->conv_sep(v2_0, _vSize, _hSize, pb[0], psb[0]);
			break;
		case 1:
			kern->conv_sep(v2_1, _vSize, _hSize, pb[2], psb[1]); //Here psb is not allocated?
			break;
		default:
			throw runtime_error("CLevel5A::compute: invalid id");
		}
	}

#if NUM_THREADS == 4

	omp_set_num_threads(4);

	#pragma omp parallel shared(L4,V2D0,V2D1,a,pb)
	{
		int id = omp_get_thread_num();

		switch(id)
		{
		case 0:
			{
				update(a[0],L4->d0s1(),V2D0,pb[0]);
				update(a[1],L4->d1s1(),V2D0,pb[0]);
				update(a[2],L4->d2s1(),V2D0,pb[0]);
				update(a[3],L4->d3s1(),V2D0,pb[0]);
			}
			break;
		case 1:
			{
				update(a[4],L4->d4s1(),V2D0,pb[0]);
				update(a[5],L4->d5s1(),V2D0,pb[0]);
				update(a[6],L4->d6s1(),V2D0,pb[0]);
				update(a[7],L4->d7s1(),V2D0,pb[0]);
			}
			break;
		case 2:
			{
				update(a[8],L4->d0s2(),V2D1,pb[2]);
				update(a[9],L4->d1s2(),V2D1,pb[2]);
				update(a[10],L4->d2s2(),V2D1,pb[2]);
				update(a[11],L4->d3s2(),V2D1,pb[2]);
			}
			break;
		case 3:
			{
				update(a[12],L4->d4s2(),V2D1,pb[2]);
				update(a[13],L4->d5s2(),V2D1,pb[2]);
				update(a[14],L4->d6s2(),V2D1,pb[2]);
				update(a[15],L4->d7s2(),V2D1,pb[2]);
			}
		}
	}

#elif NUM_THREADS == 16

	//Second do updates
	omp_set_num_threads(16);

	#pragma omp parallel shared(L4, V2D0, V2D1, a, pb)
	{
		int id = omp_get_thread_num();

		switch(id)
		{
		case 0:
			update(a[0], L4->d0s1(), V2D0, pb[0]);
			break;
		case 1:
			update(a[1], L4->d1s1(), V2D0, pb[0]);
			break;
		case 2:
			update(a[2], L4->d2s1(), V2D0, pb[0]);
			break;
		case 3:
			update(a[3],L4->d3s1(),V2D0,pb[0]);
			break;
		case 4:
			update(a[4],L4->d4s1(),V2D0,pb[0]);
			break;
		case 5:
			update(a[5],L4->d5s1(),V2D0,pb[0]);
			break;
		case 6:
			update(a[6],L4->d6s1(),V2D0,pb[0]);
			break;
		case 7:
			update(a[7],L4->d7s1(),V2D0,pb[0]);
			break;
		case 8:
			update(a[8],L4->d0s2(),V2D1,pb[2]);
			break;
		case 9:
			update(a[9],L4->d1s2(),V2D1,pb[2]);
			break;
		case 10:
			update(a[10],L4->d2s2(),V2D1,pb[2]);
			break;
		case 11:
			update(a[11],L4->d3s2(),V2D1,pb[2]);
			break;
		case 12:
			update(a[12],L4->d4s2(),V2D1,pb[2]);
			break;
		case 13:
			update(a[13],L4->d5s2(),V2D1,pb[2]);
			break;
		case 14:
			update(a[14],L4->d6s2(),V2D1,pb[2]);
			break;
		case 15:
			update(a[15],L4->d7s2(),V2D1,pb[2]);
			break;
		default:
			throw runtime_error("CLevel5A::compute: invalid direction");
		}
	}

#endif	//16 threads

#endif	//openMp
	if(keep_record())
		save();
}

double* CLevel5A::update(double* array, double* L4Input, double* V2Input, double* bf1)
{
	double* pa  = array;
	double* pL4 = L4Input;
	double* pV2 = V2Input;
	double* pb  = bf1;

	for(unsigned int i = 0; i < _hSize * _vSize; i++, pa++, pL4++, pV2++, pb++)
	{
		if(*pL4 > 0)
			*pa = *pa * _cst21 + (1 - *pa) * *pL4 * (_cst22 + _cst23 * *pV2) + _cst24 * (1 + *pa) * *pb;
		else
			*pa = *pa * _cst21 + _cst24 * (1 + *pa) * *pb;
	}
	return array;
}

void CLevel5A::get_L5A_frame()
{
	auto q = _a.begin();

	for(auto p = _input_file.begin() ; p != _input_file.end() ; p++, q++)
	{
		(*p)->read((char*) *q ,_vSize * _hSize * sizeof(double));

		if((*p)->rdstate() != ios::goodbit)
			throw ios_base::failure("CLevel5A::get_L5A_frame: unable to read frame");
	}
}

void CLevel5A::setL5ExternalInput(const char* wd,unsigned int& nbFrames, unsigned int& nbIterations,vector<string> fname)
{
	if(fname.empty())
	{
		//Use default name
		fname = vector<string>(16);
		fname[0] = "motionLevel5AD0S1.dat";	fname[1] = "motionLevel5AD1S1.dat";
		fname[2] = "motionLevel5AD2S1.dat";	fname[3] = "motionLevel5AD3S1.dat";
		fname[4] = "motionLevel5AD4S1.dat";	fname[5] = "motionLevel5AD5S1.dat";
		fname[6] = "motionLevel5AD6S1.dat";	fname[7] = "motionLevel5AD7S1.dat";

		fname[8] = "motionLevel5AD0S2.dat";	fname[9] = "motionLevel5AD1S2.dat";
		fname[10] = "motionLevel5AD2S2.dat";fname[11] = "motionLevel5AD3S2.dat";
		fname[12] = "motionLevel5AD4S2.dat";fname[13] = "motionLevel5AD5S2.dat";
		fname[14] = "motionLevel5AD6S2.dat";fname[15] = "motionLevel5AD7S2.dat";
	}

	auto q = fname.begin();
	string temp;

	_input_file = vector<ifstream*>(16);

	for(auto p = _input_file.begin() ; p != _input_file.end() ; p++, q++)
	{
		//Create full filename: assuming input files are in the working directory
		temp = wd;
		temp += "\\";
		temp += *q;

		*p = new ifstream(temp.c_str(),ios::binary | ios::in);

		(*p)->read((char*)&nbFrames,sizeof(unsigned int));
		(*p)->read((char*)&nbIterations,sizeof(unsigned int));

		unsigned int buffer;	//TODO: use seekg instead
		(*p)->read((char*)&buffer,sizeof(unsigned int));
		(*p)->read((char*)&buffer,sizeof(unsigned int));
	}
}

void CLevel5A::print_kernels(const string& output_folder, bool recreate,string optName)
{
	string base = output_folder + "/" + _level_name;
	string endName = "IsotropicKernel";
	base += endName;
	base += optName;
	string endName2 = ".txt";
	base += endName2;
	_d0s1K.print_kernel(base, recreate);
}
}
