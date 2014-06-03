#include "visio/Level4.h"
#include "visio/Level3.h"
#include <stdexcept>
#include <omp.h>

using namespace std;

namespace visio
{
CLevel4::CLevel4(int iter_per_frame, double dt,
		double A6, double C6, double D6,
		unsigned int verticalSize, unsigned int horizontalSize,
		unsigned int totFrames,
		unsigned int nbDirections,
		double sigma1J, double sigma2J,
		double sigma1K,
		double J, double K, const std::string& output_folder,
		string name, bool record, Buffer* buffer, bool testKernels)
: KernelLayer(verticalSize, horizontalSize, totFrames, nbDirections, sigma1J, sigma2J, sigma1K,
		sigma1K, name, output_folder, iter_per_frame, dt, record, buffer)
{
	_cst14 = 1-dt*A6;
	_cst15 = -dt*A6;
	_cst16 = -dt*A6*C6;
	_cst17 = dt*A6;
	_cst18 = -dt*A6*0.1*C6;
	_cst19 = -dt*A6*0.1*D6;
	_cst20 = -dt*A6*D6;

	//Delete all isotropic gaussians except for first one (m_d0s2K)
	_d1s2K.remove();
	_d2s2K.remove();
	_d3s2K.remove();

	//Weight kernels
	_d0s2K.weight(K / (2*PI*sigma1K*sigma1K) );	//Isotropic inhibitory kernel
	_d0s2K.set_intrinsic_weight(K / (2*PI*sigma1K*sigma1K) );

	_d0s1K.weight(J / (2*PI*sigma1J*sigma2J) );	//Anistropic excitatory kernels
	_d0s1K.set_intrinsic_weight(J / (2*PI*sigma1J*sigma2J) );
	_d1s1K.weight(J / (2*PI*sigma1J*sigma2J) );
	_d1s1K.set_intrinsic_weight(J / (2*PI*sigma1J*sigma2J) );
	_d2s1K.weight(J / (2*PI*sigma1J*sigma2J) );
	_d2s1K.set_intrinsic_weight(J / (2*PI*sigma1J*sigma2J) );
	_d3s1K.weight(J / (2*PI*sigma1J*sigma2J) );
	_d3s1K.set_intrinsic_weight(J / (2*PI*sigma1J*sigma2J) );

	//SVD for horizontal and vertical kernels
	_d0s1K.svd();
	_d2s1K.svd();
	_d0s2K.svd();

	//DFT for diagonal kernels
	_d1s1K.set_dft(verticalSize + _d1s1K.vsize(), horizontalSize + _d1s1K.hsize());
	_d3s1K.set_dft(verticalSize + _d3s1K.vsize(), horizontalSize + _d3s1K.hsize());

	//Create additional buffer;
#if NUM_THREADS == 4 || NUM_THREADS == 16
	_d5s1K.weight(J / (2*PI*sigma1J*sigma2J) );	//Anistropic excitatory kernels
	_d5s1K.set_intrinsic_weight(J / (2*PI*sigma1J*sigma2J) );
	_d7s1K.weight(J / (2*PI*sigma1J*sigma2J) );	//Anistropic excitatory kernels
	_d7s1K.set_intrinsic_weight(J / (2*PI*sigma1J*sigma2J) );

	_d5s1K.set_dft(verticalSize + _d1s1K.vsize(), horizontalSize + _d1s1K.hsize());
	_d7s1K.set_dft(verticalSize + _d3s1K.vsize(), horizontalSize + _d3s1K.hsize());

#if NUM_THREADS == 4 || NUM_THREADS == 16
	_d5s2K.remove();
	_d7s2K.remove();
#endif

#if NUM_THREADS == 16

	_d1s2K = _d1s1K;
	_d3s2K = _d3s1K;
	_d5s2K = _d5s1K;
	_d7s2K = _d7s1K;

	_d1s2K.set_dft(verticalSize + _d1s2K.vsize(), horizontalSize + _d1s2K.hsize());
	_d3s2K.set_dft(verticalSize + _d3s2K.vsize(), horizontalSize + _d3s2K.hsize());
	_d5s2K.set_dft(verticalSize + _d5s2K.vsize(), horizontalSize + _d5s2K.hsize());
	_d7s2K.set_dft(verticalSize + _d7s2K.vsize(), horizontalSize + _d7s2K.hsize());

#endif

#endif

	setBaseB();
	setBaseSB();
	setBaseB2();
	setBasePB();

	if(testKernels)
		time_kernel(verticalSize, horizontalSize);
}

void CLevel4::setBaseB2()
{
#if NUM_THREADS == 1
	_baseB2.push_back(_buf->_buffer2);
#elif NUM_THREADS == 4
	_baseB2.push_back(_buf->_buffer8);
	_baseB2.push_back(_buf->_buffer9);
	_baseB2.push_back(_buf->_buffer10);
	_baseB2.push_back(_buf->_buffer11);
#elif NUM_THREADS == 16
	_baseB2.push_back(_buf->_buffer32);
	_baseB2.push_back(_buf->_buffer33);
	_baseB2.push_back(_buf->_buffer34);
	_baseB2.push_back(_buf->_buffer35);
	_baseB2.push_back(_buf->_buffer36);
	_baseB2.push_back(_buf->_buffer37);
	_baseB2.push_back(_buf->_buffer38);
	_baseB2.push_back(_buf->_buffer39);
	_baseB2.push_back(_buf->_buffer40);
	_baseB2.push_back(_buf->_buffer41);
	_baseB2.push_back(_buf->_buffer42);
	_baseB2.push_back(_buf->_buffer43);
	_baseB2.push_back(_buf->_buffer44);
	_baseB2.push_back(_buf->_buffer45);
	_baseB2.push_back(_buf->_buffer46);
	_baseB2.push_back(_buf->_buffer47);
#endif
}

void CLevel4::setBaseB()
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
	_baseB.push_back(_buf->_buffer4);
	_baseB.push_back(_buf->_buffer5);
	_baseB.push_back(_buf->_buffer6);
	_baseB.push_back(_buf->_buffer7);
	_baseB.push_back(_buf->_buffer8);
	_baseB.push_back(_buf->_buffer9);
	_baseB.push_back(_buf->_buffer10);
	_baseB.push_back(_buf->_buffer11);
	_baseB.push_back(_buf->_buffer12);
	_baseB.push_back(_buf->_buffer13);
	_baseB.push_back(_buf->_buffer14);
	_baseB.push_back(_buf->_buffer15);
#endif
}

void CLevel4::setBaseSB()
{
#if NUM_THREADS == 1
	_baseSB.push_back(_buf->_buffer1);
#elif NUM_THREADS == 4
	_baseSB.push_back(_buf->_buffer4);
	_baseSB.push_back(_buf->_buffer5);
	_baseSB.push_back(_buf->_buffer6);
	_baseSB.push_back(_buf->_buffer7);
#elif NUM_THREADS == 16
	_baseSB.push_back(_buf->_buffer16);
	_baseSB.push_back(_buf->_buffer17);
	_baseSB.push_back(_buf->_buffer18);
	_baseSB.push_back(_buf->_buffer19);
	_baseSB.push_back(_buf->_buffer20);
	_baseSB.push_back(_buf->_buffer21);
	_baseSB.push_back(_buf->_buffer22);
	_baseSB.push_back(_buf->_buffer23);
	_baseSB.push_back(_buf->_buffer24);
	_baseSB.push_back(_buf->_buffer25);
	_baseSB.push_back(_buf->_buffer26);
	_baseSB.push_back(_buf->_buffer27);
	_baseSB.push_back(_buf->_buffer28);
	_baseSB.push_back(_buf->_buffer29);
	_baseSB.push_back(_buf->_buffer30);
	_baseSB.push_back(_buf->_buffer31);
#endif
}

void CLevel4::compute(void* object)
{
	CLevel3* L3 = static_cast<CLevel3*>(object);

#if NUM_THREADS ==1
	//Scale 1
	_d0s1K.conv_sep(L3->array("D0S1"),L3->vSize(),L3->hSize(),_baseB[0],_baseSB[0]);	//Anisotropic filter (FJ)
	_d0s2K.conv_sep(L3->array("D0S1"),L3->vSize(),L3->hSize(),_baseB2[0],_baseSB[0]);	//Isotropic filter (FK)
	update(m_d0s1, _baseB[0], _baseB2[0], L3->array("D4S1"),0);

	_d1s1K.conv_dft(L3->array("D1S1"),_baseB[0],L3->vSize(),L3->hSize());	//Anisotropic filter (FJ)
	_d0s2K.conv_sep(L3->array("D1S1"),L3->vSize(),L3->hSize(),_baseB2[0],_baseSB[0]);	//Isotropic filter (FK)
	update(m_d1s1, _baseB[0], _baseB2[0], L3->array("D5S1"),1);

	_d2s1K.conv_sep(L3->array("D2S1"),L3->vSize(),L3->hSize(),_baseB[0],_baseSB[0]);	//Anisotropic filter (FJ)
	_d0s2K.conv_sep(L3->array("D2S1"),L3->vSize(),L3->hSize(),_baseB2[0],_baseSB[0]);	//Isotropic filter (FK)
	update(m_d2s1, _baseB[0], _baseB2[0], L3->array("D6S1"),2);

	_d3s1K.conv_dft(L3->array("D3S1"),_baseB[0],L3->vSize(),L3->hSize());
	_d0s2K.conv_sep(L3->array("D3S1"),L3->vSize(),L3->hSize(),_baseB2[0],_baseSB[0]);	//Isotropic filter (FK)
	update(m_d3s1, _baseB[0], _baseB2[0], L3->array("D7S1"),3);

	_d0s1K.conv_sep(L3->array("D4S1"),L3->vSize(),L3->hSize(),_baseB[0],_baseSB[0]);	//Anisotropic filter (FJ)
	_d0s2K.conv_sep(L3->array("D4S1"),L3->vSize(),L3->hSize(),_baseB2[0],_baseSB[0]);	//Isotropic filter (FK)
	update(m_d4s1, _baseB[0], _baseB2[0], L3->array("D0S1"),4);

	_d1s1K.conv_dft(L3->array("D5S1"),_baseB[0],L3->vSize(),L3->hSize());
	_d0s2K.conv_sep(L3->array("D5S1"),L3->vSize(),L3->hSize(),_baseB2[0],_baseSB[0]);	//Isotropic filter (FK)
	update(m_d5s1, _baseB[0], _baseB2[0], L3->array("D1S1"),5);

	_d2s1K.conv_sep(L3->array("D6S1"),L3->vSize(),L3->hSize(),_baseB[0],_baseSB[0]);	//Anisotropic filter (FJ)
	_d0s2K.conv_sep(L3->array("D6S1"),L3->vSize(),L3->hSize(),_baseB2[0],_baseSB[0]);	//Isotropic filter (FK)
	update(m_d6s1, _baseB[0], _baseB2[0], L3->array("D2S1"),6);

	_d3s1K.conv_dft(L3->array("D7S1"),_baseB[0],L3->vSize(),L3->hSize());
	_d0s2K.conv_sep(L3->array("D7S1"),L3->vSize(),L3->hSize(),_baseB2[0],_baseSB[0]);	//Isotropic filter (FK)
	update(m_d7s1, _baseB[0], _baseB2[0], L3->array("D3S1"),7);

	//Scale 2
	_d0s1K.conv_sep(L3->array("D0S2"),L3->vSize(),L3->hSize(),_baseB[0],_baseSB[0]);	//Anisotropic filter (FJ)
	_d0s2K.conv_sep(L3->array("D0S2"),L3->vSize(),L3->hSize(),_baseB2[0],_baseSB[0]);	//Isotropic filter (FK)
	update(m_d0s2, _baseB[0], _baseB2[0], L3->array("D4S2"),0);

	_d1s1K.conv_dft(L3->array("D1S2"),_baseB[0],L3->vSize(),L3->hSize());
	_d0s2K.conv_sep(L3->array("D1S2"),L3->vSize(),L3->hSize(),_baseB2[0],_baseSB[0]);	//Isotropic filter (FK)
	update(m_d1s2, _baseB[0], _baseB2[0], L3->array("D5S2"),1);

	_d2s1K.conv_sep(L3->array("D2S2"),L3->vSize(),L3->hSize(),_baseB[0],_baseSB[0]);	//Anisotropic filter (FJ)
	_d0s2K.conv_sep(L3->array("D2S2"),L3->vSize(),L3->hSize(),_baseB2[0],_baseSB[0]);	//Isotropic filter (FK)
	update(m_d2s2, _baseB[0], _baseB2[0], L3->array("D6S2"),2);

	_d3s1K.conv_dft(L3->array("D3S2"),_baseB[0],L3->vSize(),L3->hSize());
	_d0s2K.conv_sep(L3->array("D3S2"),L3->vSize(),L3->hSize(),_baseB2[0],_baseSB[0]);	//Isotropic filter (FK)
	update(m_d3s2, _baseB[0], _baseB2[0], L3->array("D7S2"),3);

	_d0s1K.conv_sep(L3->array("D4S2"),L3->vSize(),L3->hSize(),_baseB[0],_baseSB[0]);	//Anisotropic filter (FJ)
	_d0s2K.conv_sep(L3->array("D4S2"),L3->vSize(),L3->hSize(),_baseB2[0],_baseSB[0]);	//Isotropic filter (FK)
	update(m_d4s2, _baseB[0], _baseB2[0], L3->array("D0S2"),4);

	_d1s1K.conv_dft(L3->array("D5S2"),_baseB[0],L3->vSize(),L3->hSize());
	_d0s2K.conv_sep(L3->array("D5S2"),L3->vSize(),L3->hSize(),_baseB2[0],_baseSB[0]);	//Isotropic filter (FK)
	update(m_d5s2, _baseB[0], _baseB2[0], L3->array("D1S2"),5);

	_d2s1K.conv_sep(L3->array("D6S2"),L3->vSize(),L3->hSize(),_baseB[0],_baseSB[0]);	//Anisotropic filter (FJ)
	_d0s2K.conv_sep(L3->array("D6S2"),L3->vSize(),L3->hSize(),_baseB2[0],_baseSB[0]);	//Isotropic filter (FK)
	update(m_d6s2, _baseB[0], _baseB2[0], L3->array("D2S2"),6);

	_d3s1K.conv_dft(L3->array("D7S2"),_baseB[0],L3->vSize(),L3->hSize());
	_d0s2K.conv_sep(L3->array("D7S2"),L3->vSize(),L3->hSize(),_baseB2[0],_baseSB[0]);	//Isotropic filter (FK)
	update(m_d7s2, _baseB[0], _baseB2[0], L3->array("D3S2"),7);

#else

	vector<double*> pb = _baseB;
	vector<double*> psb = _baseSB;
	vector<double*> ppb = _basePB;
	vector<double*> pb2 = _baseB2;

	auto k = _k;
	vector<double*> a = _a;

	Kernel<double>* o = &_d0s2K;

#if NUM_THREADS == 4

	omp_set_num_threads(4);

	#pragma omp parallel  shared(L3,pb,psb,ppb,pb2,k,a,o)
	{
		int id = omp_get_thread_num();

		double* input1;
		double* input2;
		double* input3;
		double* input4;

		double* pTb;
		double* pTsb;
		double* pTb2;

		Kernel<double>* dirA;
		Kernel<double>* dirB;

		double* paA;
		double* paB;
		double* paC;
		double* paD;

		double* haA;
		double* haB;
		double* haC;
		double* haD;

		unsigned int dir1, dir2;

		switch(id)
		{
		case 0:
			{
				input1 = L3->d0s1();
				input2 = L3->d1s1();
				input3 = L3->d0s2();
				input4 = L3->d1s2();

				pTb = pb[0];
				pTsb= psb[0];
				pTb2 = pb2[0];

				dirA = k[0];
				dirB = k[1];

				paA = a[0];
				paB = a[1];
				paC = a[8];
				paD = a[9];

				haA = L3->d4s1();
				haB = L3->d5s1();
				haC = L3->d4s2();
				haD = L3->d5s2();

				dir1 = 0;
				dir2 = 1;
			}
			break;
		case 1:
			{
				input1 = L3->d2s1();
				input2 = L3->d3s1();
				input3 = L3->d2s2();
				input4 = L3->d3s2();

				pTb = pb[1];
				pTsb= psb[1];
				pTb2 = pb2[1];

				dirA = k[2];
				dirB = k[3];

				paA = a[2];
				paB = a[3];
				paC = a[10];
				paD = a[11];

				haA = L3->d6s1();
				haB = L3->d7s1();
				haC = L3->d6s2();
				haD = L3->d7s2();

				dir1 = 2;
				dir2 = 3;
			}
			break;
		case 2:
			{
				input1 = L3->d4s1();
				input2 = L3->d5s1();
				input3 = L3->d4s2();
				input4 = L3->d5s2();

				pTb = pb[2];
				pTsb= psb[2];
				pTb2 = pb2[2];

				dirA = k[0];
				dirB = k[8];	//dirB = k[1];

				paA = a[4];
				paB = a[5];
				paC = a[12];
				paD = a[13];

				haA = L3->d0s1();
				haB = L3->d1s1();
				haC = L3->d0s2();
				haD = L3->d1s2();

				dir1 = 4;
				dir2 = 5;
			}
			break;
		case 3:
			{
				input1 = L3->d6s1();
				input2 = L3->d7s1();
				input3 = L3->d6s2();
				input4 = L3->d7s2();

				pTb = pb[3];
				pTsb= psb[3];
				pTb2 = pb2[3];

				dirA = k[2];
				dirB = k[9];	//dirB = k[3];

				paA = a[6];
				paB = a[7];
				paC = a[14];
				paD = a[15];

				haA = L3->d2s1();
				haB = L3->d3s1();
				haC = L3->d2s2();
				haD = L3->d3s2();

				dir1 = 6;
				dir2 = 7;
			}
		}

		dirA->conv_sep(input1,L3->vSize(),L3->hSize(),pTb,pTsb);
		o->conv_sep(input1,L3->vSize(),L3->hSize(),pTb2,pTsb);	//Isotropic filter
		update(paA,pTb,pTb2,haA,dir1);

		dirB->conv_dft(input2,pTb,L3->vSize(),L3->hSize());
		o->conv_sep(input2,L3->vSize(),L3->hSize(),pTb2,pTsb);	//Isotropic filter
		update(paB,pTb,pTb2,haB,dir2);

		dirA->conv_sep(input3,L3->vSize(),L3->hSize(),pTb,pTsb);
		o->conv_sep(input3,L3->vSize(),L3->hSize(),pTb2,pTsb);	//Isotropic filter
		update(paC,pTb,pTb2,haC,dir1);

		dirB->conv_dft(input4,pTb,L3->vSize(),L3->hSize());
		o->conv_sep(input4,L3->vSize(),L3->hSize(),pTb2,pTsb);	//Isotropic filter
		update(paD,pTb,pTb2,haD,dir2);

	}

#elif NUM_THREADS == 16

	omp_set_num_threads(16);

	#pragma omp parallel  shared(L3,pb,psb,ppb,pb2,k,a,o)
	{
		int id = omp_get_thread_num();

		double* input = nullptr;

		double* pTb = nullptr;
		double* pTsb = nullptr;
		double* pTb2 = nullptr;

		Kernel<double>* kern = nullptr;

		double* pa = nullptr;
		double* ha = nullptr;

		unsigned int dir;

		switch(id)
		{
		case 0:
			{
				input = L3->d0s1();

				pTb = pb[0];
				pTsb= psb[0];
				pTb2 = pb2[0];

				kern = k[0];

				pa = a[0];

				ha = L3->d4s1();

				dir = 0;
			}
			break;
		case 1:
			{
				input = L3->d0s2();

				pTb = pb[1];
				pTsb= psb[1];
				pTb2 = pb2[1];

				kern = k[0];

				pa = a[8];

				ha = L3->d4s2();

				dir = 0;
			}
			break;
		case 2:
			{
				input = L3->d2s1();

				pTb = pb[2];
				pTsb= psb[2];
				pTb2 = pb2[2];

				kern = k[2];

				pa = a[2];

				ha = L3->d6s1();

				dir = 2;
			}
			break;
		case 3:
			{
				input = L3->d2s2();

				pTb = pb[3];
				pTsb= psb[3];
				pTb2 = pb2[3];

				kern = k[2];

				pa = a[10];

				ha = L3->d6s2();

				dir = 2;
			}
			break;
		case 4:
			{
				input = L3->d4s1();

				pTb = pb[4];
				pTsb= psb[4];
				pTb2 = pb2[4];

				kern = k[0];

				pa = a[4];

				ha = L3->d0s1();

				dir = 4;
			}
			break;
		case 5:
			{
				input = L3->d4s2();

				pTb = pb[5];
				pTsb= psb[5];
				pTb2 = pb2[5];

				kern = k[0];

				pa = a[12];

				ha = L3->d0s2();

				dir = 4;
			}
			break;
		case 6:
			{
				input = L3->d6s1();

				pTb = pb[6];
				pTsb= psb[6];
				pTb2 = pb2[6];

				kern = k[2];

				pa = a[6];

				ha = L3->d2s1();

				dir = 6;
			}
			break;
		case 7:
			{
				input = L3->d6s2();

				pTb = pb[7];
				pTsb= psb[7];
				pTb2 = pb2[7];

				kern = k[2];

				pa = a[14];

				ha = L3->d2s2();

				dir = 6;
			}
			break;
		case 8:
			{
				input = L3->d1s1();

				pTb = pb[8];
				pTsb= psb[8];
				pTb2 = pb2[8];

				kern = k[1];

				pa = a[1];

				ha = L3->d5s1();

				dir = 1;
			}
			break;
		case 9:
			{

				input = L3->d1s2();

				pTb = pb[9];
				pTsb= psb[9];
				pTb2 = pb2[9];

				kern = k[5];

				pa = a[9];
				ha = L3->d5s2();

				dir = 1;
			}
			break;
		case 10:
			{
				input = L3->d3s1();

				pTb = pb[10];
				pTsb= psb[10];
				pTb2 = pb2[10];

				kern = k[3];

				pa = a[3];

				ha = L3->d7s1();

				dir = 3;
			}
			break;
		case 11:
			{

				input = L3->d3s2();

				pTb = pb[11];
				pTsb= psb[11];
				pTb2 = pb2[11];

				kern = k[7];

				pa = a[11];

				ha = L3->d7s2();

				dir = 3;
			}
			break;
		case 12:
			{
				input = L3->d5s1();

				pTb = pb[12];
				pTsb= psb[12];
				pTb2 = pb2[12];

				kern = k[8];	//dirB = k[1];

				pa = a[5];

				ha = L3->d1s1();

				dir = 5;
			}
			break;
		case 13:
			{

				input = L3->d5s2();

				pTb = pb[13];
				pTsb= psb[13];
				pTb2 = pb2[13];

				kern = k[10];	//dirB = k[1];

				pa = a[13];

				ha = L3->d1s2();

				dir = 5;
			}
			break;
		case 14:
			{
				input = L3->d7s1();

				pTb = pb[14];
				pTsb= psb[14];
				pTb2 = pb2[14];

				kern = k[9];	//dirB = k[3];

				pa = a[7];

				ha = L3->d3s1();

				dir = 7;
			}
			break;
		case 15:
			{

				input = L3->d7s2();

				pTb = pb[15];
				pTsb= psb[15];
				pTb2 = pb2[15];

				kern = k[11];	//dirB = k[3];

				pa = a[15];

				ha = L3->d3s2();

				dir = 7;
			}
			break;
		default:
			throw runtime_error("CLevel4::compute: invalid direction");
		}

		switch(id)
		{
		case 0:
		case 1:
		case 2:
		case 3:
		case 4:
		case 5:
		case 6:
		case 7:
			{
				kern->conv_sep(input,L3->vSize(),L3->hSize(),pTb,pTsb);
				o->conv_sep(input,L3->vSize(),L3->hSize(),pTb2,pTsb);	//Isotropic filter
				update(pa, pTb, pTb2, ha, dir);
			}
			break;
		case 8:
		case 9:
		case 10:
		case 11:
		case 12:
		case 13:
		case 14:
		case 15:
			{
				kern->conv_dft(input,pTb,L3->vSize(),L3->hSize());
				o->conv_sep(input,L3->vSize(),L3->hSize(),pTb2,pTsb);	//Isotropic filter
				update(pa,pTb,pTb2,ha,dir);
			}
			break;
		default:
			throw runtime_error("CLevel4::compute: invalid direction");
		}
	}

#endif	//16 threads

#endif	//openmp

	if(keep_record())
		save();
}

double* CLevel4::update(double* array, double* bf1, double* bf2, double* op, unsigned int direction)
{
	double *pb2 = bf2;
	double *pb1 = bf1;
	double* pa = array;
	unsigned int i,nb;
	double* pop = op;

	switch(direction)
	{
	case 0:
		{
			//First location undefined, so skip it
			pa  += 1;
			pb1 += 1;
			pop += 1;
			nb   = _vSize*_hSize - 1;
		}
		break;
	case 1:
		{
			pa  += 1;
			pb1 += 1;
			pop += 1;
			pb2 += _hSize;
			nb   = _vSize*_hSize - _hSize - 1;
		}
		break;
	case 2:
		{
			pb2 += _hSize;
			nb   = _vSize*_hSize - _hSize;
		}
		break;
	case 3:
		{
			pb2 += _hSize + 1;
			nb   = _vSize*_hSize - _hSize;
		}
		break;
	case 4:
		{
			pb2 += 1;
			nb   = _vSize*_hSize - 1;
		}
		break;
	case 5:
		{
			pa  += _hSize;
			pb1 += _hSize;
			pop += _hSize;
			pb2 += 1;
			nb   = _vSize*_hSize - _hSize -1;
		}
		break;
	case 6:
		{
			pa  += _hSize;
			pb1 += _hSize;
			pop += _hSize;
			nb = _vSize*_hSize - _hSize;
		}
		break;
	case 7:
		{
			pa  += _hSize +1;
			pb1 += _hSize +1;
			pop += _hSize +1;
			nb = _vSize*_hSize - _hSize -1;
		}
		break;
	default:
		throw runtime_error("CLevel4::update: invalid direction");
	}

	for(i = 0; i < nb; i++, pa++, pb1++, pb2++, pop++)
		*pa = _cst14 * *pa + ( _cst15 * *pb1 + _cst16 * *pb2) * *pa + _cst17 * *pb1 + _cst18 * *pb2 + _cst19 * *pop + _cst20 * *pa * *pop;

	return array;
}

void CLevel4::print_kernels(const string& output_folder, bool recreate, string optName)
{
	//1. Anisotropic kernels
	vector< Kernel<double>* > ks(8);
	ks[0] = &_d0s1K;ks[1] = &_d1s1K;ks[2]  = &_d2s1K;ks[3] = &_d3s1K;ks[4] = &_d0s1K;ks[5] = &_d1s1K;ks[6] = &_d2s1K;ks[7] = &_d3s1K;

	vector<string> ns(8);
	ns[0] = "D0S1.txt";ns[1] = "D1S1.txt";ns[2] = "D2S1.txt";ns[3] = "D3S1.txt";ns[4] = "D4S1.txt";ns[5] = "D5S1.txt";ns[6] = "D6S1.txt";ns[7] = "D7S1.txt";

	auto n = ns.begin();

	for(auto p = ks.begin() ; p != ks.end() ; p++,n++)
	{
		string base = output_folder + "/" + _level_name;
		string endName = "kernel";
		base += optName;	//Here will add an extra name tag only if optName is not default value ""
		base += endName;
		base += *n;
		(*p)->print_kernel(base, recreate);
	}

	//2. Isotropic kernel
	string base = _level_name;
	string endName = "IsotropicKernel";
	base += endName;
	base += optName;
	string endName2 = ".txt";
	base += endName2;

	_d0s2K.print_kernel(base, recreate);
}
}
