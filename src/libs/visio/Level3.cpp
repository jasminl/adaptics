#include "visio/Level3.h"
#include <stdexcept>
#include <omp.h>

using namespace std;

namespace visio
{
CLevel3::CLevel3(int iter_per_frame, double dt, double A5, double G, unsigned int verticalSize, unsigned int horizontalSize,
		unsigned int totFrames, unsigned int nbDirections, double sigma1Scale1, double sigma2Scale1,
		double sigma1Scale2, double sigma2Scale2, double theta1, double theta2,
		CLevel2& lev2, const string& output_folder, string name, bool record,
		Buffer* buffer, bool testKernels)
: KernelLayer(verticalSize, horizontalSize, totFrames, nbDirections,
		sigma1Scale1, sigma2Scale1, sigma1Scale2, sigma2Scale2,
		name, output_folder, iter_per_frame, dt, record, buffer), _theta1(theta1), _theta2(theta2)
{
	_cst13 = 1-A5*dt;

	weight_kernels(dt*A5*G);
	set_intrinsic_weight(G);

	//SVD for horizontal and vertical kernels
	_d0s1K.svd();
	_d2s1K.svd();
	_d0s2K.svd();
	_d2s2K.svd();

	//Prepare DFT for non-separable kernels
	_d1s1K.set_dft(verticalSize + _d1s1K.vsize(),horizontalSize + _d1s1K.hsize());
	_d3s1K.set_dft(verticalSize + _d3s1K.vsize(),horizontalSize + _d3s1K.hsize());
	_d1s2K.set_dft(verticalSize + _d1s2K.vsize(),horizontalSize + _d1s2K.hsize());
	_d3s2K.set_dft(verticalSize + _d3s2K.vsize(),horizontalSize + _d3s2K.hsize());

	setBaseB();
	setBaseSB();
	setBasePB();

#if NUM_THREADS == 1
	m_ks1.push_back(&_d0s1K);
	m_ks1.push_back(&_d1s1K);
	m_ks1.push_back(&_d2s1K);
	m_ks1.push_back(&_d3s1K);
	m_ks1.push_back(&_d0s1K);
	m_ks1.push_back(&_d1s1K);
	m_ks1.push_back(&_d2s1K);
	m_ks1.push_back(&_d3s1K);

	m_ks2.push_back(&_d0s2K);
	m_ks2.push_back(&_d1s2K);
	m_ks2.push_back(&_d2s2K);
	m_ks2.push_back(&_d3s2K);
	m_ks2.push_back(&_d0s2K);
	m_ks2.push_back(&_d1s2K);
	m_ks2.push_back(&_d2s2K);
	m_ks2.push_back(&_d3s2K);

	m_ns.push_back(string("E0"));
	m_ns.push_back(string("E1"));
	m_ns.push_back(string("E2"));
	m_ns.push_back(string("E3"));
	m_ns.push_back(string("E4"));
	m_ns.push_back(string("E5"));
	m_ns.push_back(string("E6"));
	m_ns.push_back(string("E7"));

	m_vs1.push_back(m_d0s1);
	m_vs1.push_back(m_d1s1);
	m_vs1.push_back(m_d2s1);
	m_vs1.push_back(m_d3s1);
	m_vs1.push_back(m_d4s1);
	m_vs1.push_back(m_d5s1);
	m_vs1.push_back(m_d6s1);
	m_vs1.push_back(m_d7s1);

	m_vs2.push_back(m_d0s2);
	m_vs2.push_back(m_d1s2);
	m_vs2.push_back(m_d2s2);
	m_vs2.push_back(m_d3s2);
	m_vs2.push_back(m_d4s2);
	m_vs2.push_back(m_d5s2);
	m_vs2.push_back(m_d6s2);
	m_vs2.push_back(m_d7s2);
#elif NUM_THREADS == 4 || NUM_THREADS == 16

		weight_extern_kernel(_d5s1K,dt*A5*G);
		weight_extern_kernel(_d7s1K,dt*A5*G);
		weight_extern_kernel(_d5s2K,dt*A5*G);
		weight_extern_kernel(_d7s2K,dt*A5*G);

		set_intrinsic_weight_to(_d5s1K,G);
		set_intrinsic_weight_to(_d7s1K,G);
		set_intrinsic_weight_to(_d5s2K,G);
		set_intrinsic_weight_to(_d7s2K,G);

		_d5s1K.set_dft(verticalSize + _d5s1K.vsize(),horizontalSize + _d5s1K.hsize());
		_d7s1K.set_dft(verticalSize + _d7s1K.vsize(),horizontalSize + _d7s1K.hsize());
		_d5s2K.set_dft(verticalSize + _d5s2K.vsize(),horizontalSize + _d5s2K.hsize());
		_d7s2K.set_dft(verticalSize + _d7s2K.vsize(),horizontalSize + _d7s2K.hsize());
#endif

	if(testKernels)
		time_kernel(verticalSize, horizontalSize);
}

void CLevel3::setBaseB()
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

void CLevel3::setBaseSB()
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
#endif
}

void CLevel3::setBasePB()
{
#if NUM_THREADS == 1
	_basePB.push_back(_buf->_buffer2);
#elif NUM_THREADS == 4
	_basePB.push_back(_buf->_buffer8);
	_basePB.push_back(_buf->_buffer9);
	_basePB.push_back(_buf->_buffer10);
	_basePB.push_back(_buf->_buffer11);
#elif NUM_THREADS == 16
	_basePB.push_back(_buf->_buffer24);
	_basePB.push_back(_buf->_buffer25);
	_basePB.push_back(_buf->_buffer26);
	_basePB.push_back(_buf->_buffer27);
	_basePB.push_back(_buf->_buffer28);
	_basePB.push_back(_buf->_buffer29);
	_basePB.push_back(_buf->_buffer30);
	_basePB.push_back(_buf->_buffer31);
#endif
}

void CLevel3::compute(void* object)
{
	CLevel2* L2 = static_cast<CLevel2*>(object);

#if NUM_THREADS==1

	vector< Kernel<double>* >::iterator p=m_ks1.begin();
	vector< Kernel<double>* >::iterator p2=m_ks2.begin();

	vector<string>::iterator q = m_ns.begin();
	vector< double* >::iterator r = m_vs1.begin();
	vector< double* >::iterator r2 = m_vs2.begin();

	unsigned int i;

	for(;p!=m_ks1.end();p++,q++,r++,p2++,r2++)
	{
		double* cv = *r;
		double* cv2= *r2;

		double* pbuf = _baseB[0];
		//1. Process even directions (0, 2, 4, 6) since they have separable kernels
		rectify(L2->array((*q).c_str()),_basePB[0],L2->vSize(),L2->hSize());
		(*p)->conv_sep(_basePB[0],L2->vSize(),L2->hSize(),_baseB[0],_baseSB[0]);

		//Scale 1
		for(i=0;i<_hSize*_vSize;i++,cv++,pbuf++)
		{
			*cv = *cv * _cst13 + *pbuf;

			//Thresholding
			if ( (*cv-_theta1) <0)	*cv = 0;
		}

		//Scale 2
		(*p2)->conv_sep(_basePB[0],L2->vSize(),L2->hSize(),_baseB[0],_baseSB[0]);

		pbuf = _baseB[0];
		for(i=0;i<_hSize*_vSize;i++,cv2++,pbuf++)
		{
			*cv2 = *cv2 * _cst13 + *pbuf;

			//Thresholding
			if ( (*cv2-_theta2) <0)	*cv2 = 0;
		}


		p++;
		p2++;
		r++;
		r2++;
		q++;

		cv = *r;
		cv2 = *r2;
		pbuf = _baseB[0];

		//2. Process odd directions (1, 3, 5, 7), since they have non-separable kernels
		rectify(L2->array((*q).c_str()),_basePB[0],L2->vSize(),L2->hSize());
		(*p)->conv_dft(_basePB[0],_baseB[0],L2->vSize(),L2->hSize());

		for(i=0;i<_hSize*_vSize;i++,cv++,pbuf++)
		{
			*cv = *cv * _cst13 + *pbuf;

			//Thresholding
			if ( (*cv-_theta1) <0)	*cv = 0;
		}

		pbuf = _baseB[0];

		(*p2)->conv_dft(_basePB[0],_baseB[0],L2->vSize(),L2->hSize());


		for(i=0;i<_hSize*_vSize;i++,cv2++,pbuf++)
		{
			*cv2 = *cv2 * _cst13 + *pbuf;

			//Thresholding
			if ( (*cv2-_theta2) <0)	*cv2 = 0;
		}
	}
#else

	double* cv = nullptr;
	double* pbuf = nullptr;
	unsigned int i;

	double theta1 = _theta1;

	vector<double*> pb = _baseB;
	vector<double*> psb = _baseSB;
	vector<double*> ppb = _basePB;

	auto k = _k;
	vector<double*> a = _a;

	double cst13 = _cst13;
	unsigned int size = _hSize * _vSize;

#if NUM_THREADS == 4
	double theta2 = _theta2;

	omp_set_num_threads(4);

	#pragma omp parallel private(i,cv,pbuf) shared(L2,pb,psb,ppb,k,a,theta1,theta2,cst13,size)
	{
		int id = omp_get_thread_num();

		double* input1;
		double* input2;
		double* pTb;
		double* pTsb;
		double* pTpb;
		Kernel<double>* dirAS1;
		Kernel<double>* dirBS1;
		Kernel<double>* dirAS2;
		Kernel<double>* dirBS2;

		double* paA;
		double* paB;
		double* paC;
		double* paD;

		switch(id)
		{
		case 0:
			{
				input1 = L2->e0();
				input2 = L2->e1();

				pTb = pb[0];
				pTsb= psb[0];
				pTpb= ppb[0];

				dirAS1 = k[0];
				dirBS1 = k[1];
				dirAS2 = k[4];
				dirBS2 = k[5];

				paA = a[0];
				paB = a[8];
				paC = a[1];
				paD = a[9];
			}
			break;
		case 1:
			{
				input1 = L2->e2();
				input2 = L2->e3();

				pTb = pb[1];
				pTsb= psb[1];
				pTpb= ppb[1];

				dirAS1 = k[2];
				dirBS1 = k[3];
				dirAS2 = k[6];
				dirBS2 = k[7];

				paA = a[2];
				paB = a[10];
				paC = a[3];
				paD = a[11];
			}
			break;
		case 2:
			{
				input1 = L2->e4();
				input2 = L2->e5();

				pTb = pb[2];
				pTsb= psb[2];
				pTpb= ppb[2];

				dirAS1 = k[0];
				dirBS1 = k[8];	//dirBS1 = k[1];
				dirAS2 = k[4];
				dirBS2 = k[10];	//dirBS2 = k[5];

				paA = a[4];
				paB = a[12];
				paC = a[5];
				paD = a[13];
			}
			break;
		case 3:
			{
				input1 = L2->e6();
				input2 = L2->e7();

				pTb = pb[3];
				pTsb= psb[3];
				pTpb= ppb[3];

				dirAS1 = k[2];
				dirBS1 = k[9];	//dirBS1 = k[3];
				dirAS2 = k[6];
				dirBS2 = k[11];	//dirBS2 = k[7];

				paA = a[6];
				paB = a[14];
				paC = a[7];
				paD = a[15];
			}
		}

		rectify(input1,pTpb,L2->vSize(),L2->hSize());
		dirAS1->conv_sep(pTpb,L2->vSize(),L2->hSize(),pTb,pTsb);

		cv = paA;
		pbuf = pTb;
		for(i=0;i<size;i++,cv++,pbuf++)
		{
			*cv = *cv * cst13 + *pbuf;

			if ( (*cv-theta1) <0)	*cv = 0;
		}

		dirAS2->conv_sep(pTpb,L2->vSize(),L2->hSize(),pTb,pTsb);

		cv = paB;
		pbuf = pTb;
		for(i=0;i<size;i++,cv++,pbuf++)
		{
			*cv = *cv * cst13 + *pbuf;

			if ( (*cv-theta2) <0)	*cv = 0;
		}

		rectify(input2,pTpb,L2->vSize(),L2->hSize());

		dirBS1->conv_dft(pTpb,pTb,L2->vSize(),L2->hSize());

		cv = paC;
		pbuf = pTb;
		for(i=0;i<size;i++,cv++,pbuf++)
		{
			*cv = *cv * cst13 + *pbuf;

			if ( (*cv-theta1) <0)	*cv = 0;
		}

		dirBS2->conv_dft(pTpb,pTb,L2->vSize(),L2->hSize());

		cv = paD;
		pbuf = pTb;
		for(i=0;i<size;i++,cv++,pbuf++)
		{
			*cv = *cv * cst13 + *pbuf;

			if ( (*cv-theta2) <0)	*cv = 0;
		}

	}

#elif NUM_THREADS == 16
	omp_set_num_threads(16);

	#pragma omp parallel private(i, cv, pbuf) shared(L2, pb, psb, ppb, k, a, theta1/*, theta2*/, cst13, size)
	{
		int id = omp_get_thread_num();

		double* input;
		double* pTb;
		double* pTsb;
		double* pTpb;
		Kernel<double>* kern = nullptr;

		unsigned int chopSize = chop_input2(size);		//Half-size of input array

		switch(id)
		{
		case 0:
			{
				input = L2->e0();

				pTb = pb[0];
				pTsb= psb[0];
				pTpb= ppb[0];

				kern = k[0];

				cv = a[0];

				rectify(input, pTpb, chopSize);
			}
			break;
		case 1:
			{
				input = L2->e0();

				pTb = pb[8];
				pTsb= psb[1];
				pTpb= ppb[0];

				kern = k[4];

				cv = a[8];

				rectify(input + chopSize, pTpb + chopSize, size - chopSize);
			}
			break;
		case 2:
			{
				input = L2->e2();

				pTb = pb[2];
				pTsb= psb[2];
				pTpb= ppb[2];

				kern = k[2];

				cv = a[2];

				rectify(input, pTpb, size);
			}
			break;
		case 3:
			{
				input = L2->e2();

				pTb = pb[10];
				pTsb= psb[3];
				pTpb= ppb[2];

				kern = k[6];

				cv = a[10];

				rectify(input + chopSize, pTpb + chopSize, size - chopSize);
			}
			break;
		case 4:
			{
				input = L2->e4();

				pTb = pb[4];
				pTsb= psb[4];
				pTpb= ppb[4];

				kern = k[0];

				cv = a[4];

				rectify(input, pTpb, size);
			}
			break;
		case 5:
			{
				input = L2->e4();

				pTb = pb[12];
				pTsb= psb[5];
				pTpb= ppb[4];

				kern = k[4];

				cv = a[12];

				rectify(input + chopSize, pTpb + chopSize, size - chopSize);
			}
			break;
		case 6:
			{
				input = L2->e6();

				pTb = pb[6];
				pTsb= psb[6];
				pTpb= ppb[6];

				kern = k[2];

				cv = a[6];

				rectify(input, pTpb, size);
			}
			break;
		case 7:
			{
				input = L2->e6();

				pTb = pb[14];
				pTsb= psb[7];
				pTpb= ppb[6];

				kern = k[6];

				cv = a[14];

				rectify(input + chopSize, pTpb + chopSize, size - chopSize);
			}
			break;
		case 8:
			{
				input = L2->e1();

				pTb = pb[1];
				pTpb= ppb[1];

				kern = k[1];

				cv = a[1];

				rectify(input, pTpb, chopSize);
			}
			break;
		case 9:
			{
				input = L2->e1();

				pTb = pb[9];
				pTpb= ppb[1];

				kern = k[5];

				cv = a[9];

				rectify(input + chopSize, pTpb + chopSize, size - chopSize);
			}
			break;
		case 10:
			{
				input = L2->e3();

				pTb = pb[3];
				pTpb= ppb[3];

				kern = k[3];

				cv = a[3];

				rectify(input, pTpb, size);
			}
			break;
		case 11:
			{
				input = L2->e3();

				pTb = pb[11];
				pTpb= ppb[3];

				kern = k[7];

				cv = a[11];

				rectify(input + chopSize, pTpb + chopSize, size - chopSize);
			}
			break;
		case 12:
			{
				input = L2->e5();

				pTb = pb[5];
				pTpb= ppb[5];

				kern = k[8];	//dirBS1 = k[1];

				cv = a[5];

				rectify(input, pTpb, size);
			}
			break;
		case 13:
			{
				input = L2->e5();

				pTb = pb[13];
				pTpb= ppb[5];

				kern = k[10];	//dirBS2 = k[5];

				cv = a[13];

				rectify(input + chopSize, pTpb + chopSize, size - chopSize);
			}
			break;
		case 14:
			{
				input = L2->e7();

				pTb = pb[7];
				pTpb= ppb[7];

				kern = k[9];	//dirBS1 = k[3];

				cv = a[7];

				rectify(input, pTpb, size);
			}
			break;
		case 15:
			{
				input = L2->e7();

				pTb = pb[15];
				pTpb= ppb[7];

				kern = k[11];	//dirBS2 = k[7];

				cv = a[15];

				rectify(input + chopSize, pTpb + chopSize, size - chopSize);
			}
			break;
		default:
			throw runtime_error("CLevel3::compute: invalid direction");
		}

#pragma omp barrier		//Wait for all arrays to have been rectified

		switch (id)
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
				kern->conv_sep(pTpb,L2->vSize(),L2->hSize(),pTb,pTsb);
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
				kern->conv_dft(pTpb,pTb,L2->vSize(),L2->hSize());
			}
			break;
		default:
			throw runtime_error("CLayer3::compute: invalid direction specified");
		}

		pbuf = pTb;

		for(i = 0; i < size; i++, cv++, pbuf++)
		{
			*cv = *cv * cst13 + *pbuf;

			if ((*cv - theta1) < 0)
				*cv = 0;
		}
	}

#endif	//16 Threads

#endif	//Multiple threads

	if(keep_record())
		save();
}
}
