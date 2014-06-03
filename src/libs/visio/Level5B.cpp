#include "visio/Level5B.h"
#include "visio/Level5A.h"
#include <stdexcept>
#include <omp.h>

using namespace std;

namespace visio
{
CLevel5B::CLevel5B(int iter_per_frame, double dt,
		double A8, double D8,
		unsigned int verticalSize, unsigned int horizontalSize,
		unsigned int totFrames,
		unsigned int nbDirections,
		double alpha,
		double sigmaxs1, double sigmays1,double sigmaPs1,double Ls1,double thetas1,
		double sigmaxs2, double sigmays2,double sigmaPs2,double Ls2,double thetas2,
		double minW, double maxW, unsigned int type, const std::string& output_folder,
		string name, bool record, Buffer* buffer, bool testKernels)
: KernelLayer(verticalSize, horizontalSize, totFrames, nbDirections,
		sigmaxs1, sigmays1, sigmaxs2, sigmays2, name,
		output_folder, iter_per_frame, dt, record, buffer), _thetas1(thetas1), _thetas2(thetas2)
{
	//Calculate constants
	_cst25 = 1-dt*A8;
	_cst26 = dt*A8;
	_cst30 = dt*A8*alpha;

	//Set CLevel6 Pointer to NULL
	_L6 = nullptr;

	//Weight filters
	_p_s1 = Kernel<double>(0,sigmaPs1,sigmaPs1,THRES);
	_p_s2 = Kernel<double>(0,sigmaPs2,sigmaPs2,THRES);
	_p_s1.weight( (-dt*A8*D8) / (2 * PI * sigmaPs1*sigmaPs1));
	_p_s1.set_intrinsic_weight(1/ (2 * PI * sigmaPs1*sigmaPs1));
	_p_s2.weight( (-dt*A8*D8) / (2 * PI * sigmaPs2*sigmaPs2));
	_p_s2.set_intrinsic_weight(1/ (2 * PI * sigmaPs2*sigmaPs2));

	_d0s1K.weight(Ls1 / (2*PI * sigmaxs1 * sigmays1));
	_d0s1K.set_intrinsic_weight(Ls1 / (2*PI * sigmaxs1 * sigmays1));
	_d1s1K.weight(Ls1 / (2*PI * sigmaxs1 * sigmays1));
	_d1s1K.set_intrinsic_weight(Ls1 / (2*PI * sigmaxs1 * sigmays1));
	_d2s1K.weight(Ls1 / (2*PI * sigmaxs1 * sigmays1));
	_d2s1K.set_intrinsic_weight(Ls1 / (2*PI * sigmaxs1 * sigmays1));
	_d3s1K.weight(Ls1 / (2*PI * sigmaxs1 * sigmays1));
	_d3s1K.set_intrinsic_weight(Ls1 / (2*PI * sigmaxs1 * sigmays1));

	_d0s2K.weight(Ls2 / (2*PI * sigmaxs2 * sigmays2));
	_d0s2K.set_intrinsic_weight(Ls2 / (2*PI * sigmaxs2 * sigmays2));
	_d1s2K.weight(Ls2 / (2*PI * sigmaxs2 * sigmays2));
	_d1s2K.set_intrinsic_weight(Ls2 / (2*PI * sigmaxs2 * sigmays2));
	_d2s2K.weight(Ls2 / (2*PI * sigmaxs2 * sigmays2));
	_d2s2K.set_intrinsic_weight(Ls2 / (2*PI * sigmaxs2 * sigmays2));
	_d3s2K.weight(Ls2 / (2*PI * sigmaxs2 * sigmays2));
	_d3s2K.set_intrinsic_weight(Ls2 / (2*PI * sigmaxs2 * sigmays2));

	//SVD for horizontal, vertical and isotropic kernels
	_d0s1K.svd();
	_d2s1K.svd();
	_d0s2K.svd();
	_d2s2K.svd();
	_p_s1.svd();
	_p_s2.svd();

	//FFT for diagonal kernels
	_d1s1K.set_dft(verticalSize + _d1s1K.vsize(),horizontalSize + _d1s1K.hsize());
	_d3s1K.set_dft(verticalSize + _d3s1K.vsize(),horizontalSize + _d3s1K.hsize());
	_d1s2K.set_dft(verticalSize + _d1s2K.vsize(),horizontalSize + _d1s2K.hsize());
	_d3s2K.set_dft(verticalSize + _d3s2K.vsize(),horizontalSize + _d3s2K.hsize());

	//Determine directional gradient for cross-direction competition w^de
    _dg = vector<double>(2 * nbDirections);
	auto p = _dg.begin();
	if (type == 0)
	{
		//Circular kernel
		vector<double> shape(nbDirections);
		shape[0] = 0; shape[1] = 1.0/4; shape[2] = 1.0/2; shape[3] = 3.0/4;
		shape[4] = 1; shape[5] = 3.0/4; shape[6] = 1.0/2; shape[7] = 1.0/4;

		for(auto q = shape.begin(); q != shape.end(); p++,q++)
			*p = *q * (maxW - minW) + minW;
	}
	else
	{
		//Step kernel
		_dg[0] = 0;
		fill(_dg.begin() + 1, _dg.begin() + nbDirections, maxW);
	}

	copy(_dg.begin(), _dg.begin() + nbDirections, _dg.begin() + nbDirections);

#if NUM_THREADS==4 || NUM_THREADS == 16

	/*Note: the extra kernels are meant to avoid race conditions when using k.conv_dft in CLevel5B::update.
	This is not necessary for separable kernels, since they don't use reserved FFTW space*/
	_d5s1K.weight(Ls1 / (2*PI * sigmaxs1 * sigmays1));
	_d5s1K.set_intrinsic_weight(Ls1 / (2*PI * sigmaxs1 * sigmays1));
	_d7s1K.weight(Ls1 / (2*PI * sigmaxs1 * sigmays1));
	_d7s1K.set_intrinsic_weight(Ls1 / (2*PI * sigmaxs1 * sigmays1));

	_d5s1K.set_dft(verticalSize + _d1s1K.vsize(),horizontalSize + _d1s1K.hsize());
	_d7s1K.set_dft(verticalSize + _d3s1K.vsize(),horizontalSize + _d3s1K.hsize());

	_d5s2K.weight(Ls2 / (2*PI * sigmaxs2 * sigmays2));
	_d5s2K.set_intrinsic_weight(Ls2 / (2*PI * sigmaxs2 * sigmays2));
	_d7s2K.weight(Ls2 / (2*PI * sigmaxs2 * sigmays2));
	_d7s2K.set_intrinsic_weight(Ls2 / (2*PI * sigmaxs2 * sigmays2));

	_d5s2K.set_dft(verticalSize + _d1s2K.vsize(),horizontalSize + _d1s2K.hsize());
	_d7s2K.set_dft(verticalSize + _d3s2K.vsize(),horizontalSize + _d3s2K.hsize());

#endif

	setBaseB();
	setBaseSB();
	setBasePB();
	setBaseV();

	_level_name = "motionLevel5B";

	if(testKernels == true)
	{
		time_kernel(verticalSize, horizontalSize);
		time_kernel(verticalSize, horizontalSize, &_p_s1, "P_scale1");
		time_kernel(verticalSize, horizontalSize, &_p_s2, "P_scale2");
	}
}

void CLevel5B::setBaseB()
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

void CLevel5B::setBaseSB()
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

void CLevel5B::setBasePB()
{
#if NUM_THREADS == 1
	_basePB.push_back(_buf->_buffer2);
#elif NUM_THREADS == 4
	_basePB.push_back(_buf->_buffer8);
	_basePB.push_back(_buf->_buffer9);
	_basePB.push_back(_buf->_buffer10);
	_basePB.push_back(_buf->_buffer11);
#elif NUM_THREADS == 16
	_basePB.push_back(_buf->_buffer32);
	_basePB.push_back(_buf->_buffer33);
	_basePB.push_back(_buf->_buffer34);
	_basePB.push_back(_buf->_buffer35);
	_basePB.push_back(_buf->_buffer36);
	_basePB.push_back(_buf->_buffer37);
	_basePB.push_back(_buf->_buffer38);
	_basePB.push_back(_buf->_buffer39);
	_basePB.push_back(_buf->_buffer40);
	_basePB.push_back(_buf->_buffer41);
	_basePB.push_back(_buf->_buffer42);
	_basePB.push_back(_buf->_buffer43);
	_basePB.push_back(_buf->_buffer44);
	_basePB.push_back(_buf->_buffer45);
	_basePB.push_back(_buf->_buffer46);
	_basePB.push_back(_buf->_buffer47);
#endif
}

void CLevel5B::setBaseV()
{
#if NUM_THREADS == 1
	_baseV.push_back(_buf->_buffer3);
	_baseV.push_back(_buf->_buffer4);
	_baseV.push_back(_buf->_buffer5);
	_baseV.push_back(_buf->_buffer6);
	_baseV.push_back(_buf->_buffer7);
	_baseV.push_back(_buf->_buffer8);
	_baseV.push_back(_buf->_buffer9);
	_baseV.push_back(_buf->_buffer10);
#elif NUM_THREADS == 4
	_baseV.push_back(_buf->_buffer12);
	_baseV.push_back(_buf->_buffer13);
	_baseV.push_back(_buf->_buffer14);
	_baseV.push_back(_buf->_buffer15);
	_baseV.push_back(_buf->_buffer16);
	_baseV.push_back(_buf->_buffer17);
	_baseV.push_back(_buf->_buffer18);
	_baseV.push_back(_buf->_buffer19);
#elif NUM_THREADS == 16
	_baseV.push_back(_buf->_buffer48);
	_baseV.push_back(_buf->_buffer49);
	_baseV.push_back(_buf->_buffer50);
	_baseV.push_back(_buf->_buffer51);
	_baseV.push_back(_buf->_buffer52);
	_baseV.push_back(_buf->_buffer53);
	_baseV.push_back(_buf->_buffer54);
	_baseV.push_back(_buf->_buffer55);
	_baseV.push_back(_buf->_buffer56);
	_baseV.push_back(_buf->_buffer57);
	_baseV.push_back(_buf->_buffer58);
	_baseV.push_back(_buf->_buffer59);
	_baseV.push_back(_buf->_buffer60);
	_baseV.push_back(_buf->_buffer61);
	_baseV.push_back(_buf->_buffer62);
	_baseV.push_back(_buf->_buffer63);
#endif
}

void CLevel5B::compute(void* object)
{
	CLevel5A* L5A = static_cast<CLevel5A*>(object);

#if NUM_THREADS==1

	//1. Convolve T once for scale 1
	_p_s1.conv_sepRectified(_L6->array("D0S1"),_vSize,_hSize,_baseV[0],_baseSB[0]);
	_p_s1.conv_sepRectified(_L6->array("D1S1"),_vSize,_hSize,_baseV[1],_baseSB[0]);
	_p_s1.conv_sepRectified(_L6->array("D2S1"),_vSize,_hSize,_baseV[2],_baseSB[0]);
	_p_s1.conv_sepRectified(_L6->array("D3S1"),_vSize,_hSize,_baseV[3],_baseSB[0]);
	_p_s1.conv_sepRectified(_L6->array("D4S1"),_vSize,_hSize,_baseV[4],_baseSB[0]);
	_p_s1.conv_sepRectified(_L6->array("D5S1"),_vSize,_hSize,_baseV[5],_baseSB[0]);
	_p_s1.conv_sepRectified(_L6->array("D6S1"),_vSize,_hSize,_baseV[6],_baseSB[0]);
	_p_s1.conv_sepRectified(_L6->array("D7S1"),_vSize,_hSize,_baseV[7],_baseSB[0]);

	//2. Update qties at scale 1
	updateSep(m_d0s1,0,L5A->array("D0S1"),_d0s1K,_thetas1,_L6->array("D0S1"));
	update(m_d1s1,1,L5A->array("D1S1"),_d1s1K,_thetas1,_L6->array("D1S1"));
	updateSep(m_d2s1,2,L5A->array("D2S1"),_d2s1K,_thetas1,_L6->array("D2S1"));
	update(m_d3s1,3,L5A->array("D3S1"),_d3s1K,_thetas1,_L6->array("D3S1"));
	updateSep(m_d4s1,4,L5A->array("D4S1"),_d0s1K,_thetas1,_L6->array("D4S1"));
	update(m_d5s1,5,L5A->array("D5S1"),_d1s1K,_thetas1,_L6->array("D5S1"));
	updateSep(m_d6s1,6,L5A->array("D6S1"),_d2s1K,_thetas1,_L6->array("D6S1"));
	update(m_d7s1,7,L5A->array("D7S1"),_d3s1K,_thetas1,_L6->array("D7S1"));

	//3. Convolve T once for scale 2
	_p_s2.conv_sepRectified(_L6->array("D0S2"),_vSize,_hSize,_baseV[0],_baseSB[0]);
	_p_s2.conv_sepRectified(_L6->array("D1S2"),_vSize,_hSize,_baseV[1],_baseSB[0]);
	_p_s2.conv_sepRectified(_L6->array("D2S2"),_vSize,_hSize,_baseV[2],_baseSB[0]);
	_p_s2.conv_sepRectified(_L6->array("D3S2"),_vSize,_hSize,_baseV[3],_baseSB[0]);
	_p_s2.conv_sepRectified(_L6->array("D4S2"),_vSize,_hSize,_baseV[4],_baseSB[0]);
	_p_s2.conv_sepRectified(_L6->array("D5S2"),_vSize,_hSize,_baseV[5],_baseSB[0]);
	_p_s2.conv_sepRectified(_L6->array("D6S2"),_vSize,_hSize,_baseV[6],_baseSB[0]);
	_p_s2.conv_sepRectified(_L6->array("D7S2"),_vSize,_hSize,_baseV[7],_baseSB[0]);

	//4. Update qties at scale 2
	updateSep(m_d0s2,0,L5A->array("D0S2"),_d0s2K,_thetas2,_L6->array("D0S2"));
	update(m_d1s2,1,L5A->array("D1S2"),_d1s2K,_thetas2,_L6->array("D1S2"));
	updateSep(m_d2s2,2,L5A->array("D2S2"),_d2s2K,_thetas2,_L6->array("D2S2"));
	update(m_d3s2,3,L5A->array("D3S2"),_d3s2K,_thetas2,_L6->array("D3S2"));
	updateSep(m_d4s2,4,L5A->array("D4S2"),_d0s2K,_thetas2,_L6->array("D4S2"));
	update(m_d5s2,5,L5A->array("D5S2"),_d1s2K,_thetas2,_L6->array("D5S2"));
	updateSep(m_d6s2,6,L5A->array("D6S2"),_d2s2K,_thetas2,_L6->array("D6S2"));
	update(m_d7s2,7,L5A->array("D7S2"),_d3s2K,_thetas2,_L6->array("D7S2"));

#else

	Kernel<double>* pOs1 = &_p_s1;
	Kernel<double>* pOs2 = &_p_s2;

	CLevel6* p6 = _L6;

	vector<double*> pb = _baseB;
	vector<double*> psb = _baseSB;
	vector<double*> ppb = _basePB;
	vector<double*> v = _baseV;

	vector< Kernel<double>* > k = _k;
	vector<double*> a = _a;

	CLevel5B* pObject = this;

#if NUM_THREADS == 4
	omp_set_num_threads(4);

	//Parrallelize cvd computation
	#pragma omp parallel  shared(L5A,pOs1,pOs2,p6,pb,ppb,psb,v,a,k,pObject)
	{
		int id = omp_get_thread_num();

		double* _ppb;
		double* _psb;
		double* _pb;

		double* _vA;
		double* _vB;
		double* input1;
		double* input2;
		double* input3;
		double* input4;

		double* L6input1;
		double* L6input2;
		double* L6input3;
		double* L6input4;

		double* q1;
		double* q2;
		double* q3;
		double* q4;

		Kernel<double>* kA;
		Kernel<double>* kB;
		Kernel<double>* kC;
		Kernel<double>* kD;

		unsigned int dirA;
		unsigned int dirB;

		switch(id)
		{
		case 0:
			{

				_ppb = ppb[0];
				_psb = psb[0];
				_pb = pb[0];

				_vA = v[0];
				_vB = v[1];

				L6input1 = p6->d0s1();
				L6input2 = p6->d1s1();
				L6input3 = p6->d0s2();
				L6input4 = p6->d1s2();

				q1 = a[0];
				q2 = a[1];
				q3 = a[8];
				q4 = a[9];

				dirA = 0;
				dirB = 1;

				input1 = L5A->d0s1();
				input2 = L5A->d1s1();
				input3 = L5A->d0s2();
				input4 = L5A->d1s2();

				kA = k[0];
				kB = k[1];
				kC = k[4];
				kD = k[5];
			}
			break;
		case 1:
			{

				_ppb = ppb[1];
				_psb = psb[1];
				_pb = pb[1];

				_vA = v[2];
				_vB = v[3];

				L6input1 = p6->d2s1();
				L6input2 = p6->d3s1();
				L6input3 = p6->d2s2();
				L6input4 = p6->d3s2();

				q1 = a[2];
				q2 = a[3];
				q3 = a[10];
				q4 = a[11];

				dirA = 2;
				dirB = 3;

				input1 = L5A->d2s1();
				input2 = L5A->d3s1();
				input3 = L5A->d2s2();
				input4 = L5A->d3s2();

				kA = k[2];
				kB = k[3];
				kC = k[6];
				kD = k[7];
			}
			break;
		case 2:
			{

				_ppb = ppb[2];
				_psb = psb[2];
				_pb = pb[2];

				_vA = v[4];
				_vB = v[5];

				L6input1 = p6->d4s1();
				L6input2 = p6->d5s1();
				L6input3 = p6->d4s2();
				L6input4 = p6->d5s2();

				q1 = a[4];
				q2 = a[5];
				q3 = a[12];
				q4 = a[13];

				dirA = 4;
				dirB = 5;

				input1 = L5A->d4s1();
				input2 = L5A->d5s1();
				input3 = L5A->d4s2();
				input4 = L5A->d5s2();

				kA = k[0];
				kB = k[8];	//kB = k[1];
				kC = k[4];
				kD = k[10]; //kD = k[5];
			}
			break;
		case 3:
			{
				_ppb = ppb[3];
				_psb = psb[3];
				_pb = pb[3];

				_vA = v[6];
				_vB = v[7];

				L6input1 = p6->d6s1();
				L6input2 = p6->d7s1();
				L6input3 = p6->d6s2();
				L6input4 = p6->d7s2();

				q1 = a[6];
				q2 = a[7];
				q3 = a[14];
				q4 = a[15];

				dirA = 6;
				dirB = 7;

				input1 = L5A->d6s1();
				input2 = L5A->d7s1();
				input3 = L5A->d6s2();
				input4 = L5A->d7s2();

				kA = k[2];
				kB = k[9];	//kB = k[3];
				kC = k[6];
				kD = k[11];	//kD = k[7];
			}
		}

		rectify(L6input1,_ppb,p6->vSize(),p6->hSize());
		pOs1->conv_sep(_ppb,p6->vSize(),p6->hSize(),_vA,_psb);

		rectify(L6input2,_ppb,p6->vSize(),p6->hSize());
		pOs1->conv_sep(_ppb,p6->vSize(),p6->hSize(),_vB,_psb);

		#pragma omp barrier

		updateSep(q1,dirA,input1,*kA,pObject,_pb,_psb,_ppb,pObject->_thetas1,L6input1);
		update(q2,dirB,input2,*kB,pObject,_pb,_ppb,pObject->_thetas1,L6input2);

		#pragma omp barrier

		rectify(L6input3,_ppb,p6->vSize(),p6->hSize());
		pOs2->conv_sep(_ppb,p6->vSize(),p6->hSize(),_vA,_psb);

		rectify(L6input4,_ppb,p6->vSize(),p6->hSize());
		pOs2->conv_sep(_ppb,p6->vSize(),p6->hSize(),_vB,_psb);

		#pragma omp barrier

		updateSep(q3,dirA,input3,*kC,pObject,_pb,_psb,_ppb,pObject->_thetas2,L6input3);
		update(q4,dirB,input4,*kD,pObject,_pb,_ppb,pObject->_thetas2,L6input4);

	}

#elif NUM_THREADS == 16

	omp_set_num_threads(16);

	//Parrallelize cvd computation
	#pragma omp parallel  shared(L5A,pOs1,pOs2,p6,pb,ppb,psb,v,a,k,pObject)
	{
		int id = omp_get_thread_num();

		double* _ppb;
		double* _psb;
		double* _pb;

		double* _v;
		double* input;

		double* L6input;

		double* q;

		Kernel<double>* kern;

		unsigned int dir;

		switch(id)
		{
		case 0:
			{
				_ppb = ppb[0];
				_psb = psb[0];
				_pb = pb[0];

				_v = v[0];

				L6input = p6->d0s1();

				q = a[0];

				dir = 0;

				input = L5A->d0s1();

				kern = k[0];
			}
			break;
		case 1:
			{
				_ppb = ppb[1];
				_psb = psb[1];
				_pb = pb[1];

				_v = v[8];

				L6input = p6->d0s2();

				q = a[8];

				dir = 0;

				input = L5A->d0s2();

				kern = k[4];
			}
			break;
		case 2:
			{

				_ppb = ppb[2];
				_psb = psb[2];
				_pb = pb[2];

				_v = v[2];

				L6input = p6->d2s1();

				q = a[2];

				dir = 2;

				input = L5A->d2s1();

				kern = k[2];
			}
			break;
		case 3:
			{

				_ppb = ppb[3];
				_psb = psb[3];
				_pb = pb[3];

				_v = v[10];

				L6input = p6->d2s2();

				q = a[10];

				dir = 2;

				input = L5A->d2s2();

				kern = k[6];
			}
			break;
		case 4:
			{

				_ppb = ppb[4];
				_psb = psb[4];
				_pb = pb[4];

				_v = v[4];

				L6input = p6->d4s1();

				q = a[4];

				dir = 4;

				input = L5A->d4s1();

				kern = k[0];
			}
			break;
		case 5:
			{

				_ppb = ppb[5];
				_psb = psb[5];
				_pb = pb[5];

				_v = v[12];

				L6input = p6->d4s2();

				q = a[12];

				dir = 4;

				input = L5A->d4s2();

				kern = k[4];
			}
			break;
		case 6:
			{
				_ppb = ppb[6];
				_psb = psb[6];
				_pb = pb[6];

				_v = v[6];

				L6input = p6->d6s1();

				q = a[6];

				dir = 6;

				input = L5A->d6s1();

				kern = k[2];

			}
			break;
		case 7:
			{
				_ppb = ppb[7];
				_psb = psb[7];
				_pb = pb[7];

				_v = v[14];

				L6input = p6->d6s2();

				q = a[14];

				dir = 6;

				input = L5A->d6s2();

				kern = k[6];
			}
			break;
		case 8:
			{
				_ppb = ppb[8];
				_psb = psb[8];
				_pb = pb[8];

				_v = v[1];

				L6input = p6->d1s1();

				q = a[1];

				dir = 1;

				input = L5A->d1s1();

				kern = k[1];

			}
			break;
		case 9:
			{
				_ppb = ppb[9];
				_psb = psb[9];
				_pb = pb[9];

				_v = v[9];

				L6input = p6->d1s2();

				q = a[9];

				dir = 1;

				input = L5A->d1s2();

				kern = k[5];
			}
			break;
		case 10:
			{

				_ppb = ppb[10];
				_psb = psb[10];
				_pb = pb[10];

				_v = v[3];

				L6input = p6->d3s1();

				q = a[3];

				dir = 3;

				input = L5A->d3s1();

				kern = k[3];

			}
			break;
		case 11:
			{

				_ppb = ppb[11];
				_psb = psb[11];
				_pb = pb[11];

				_v = v[11];

				L6input = p6->d3s2();

				q = a[11];

				dir = 3;

				input = L5A->d3s2();

				kern = k[7];
			}
			break;
		case 12:
			{

				_ppb = ppb[12];
				_psb = psb[12];
				_pb = pb[12];

				_v = v[5];

				L6input = p6->d5s1();

				q = a[5];

				dir = 5;

				input = L5A->d5s1();

				kern = k[8];	//kB = k[1];

			}
			break;
		case 13:
			{

				_ppb = ppb[13];
				_psb = psb[13];
				_pb = pb[13];

				_v = v[13];

				L6input = p6->d5s2();

				q = a[13];

				dir = 5;

				input = L5A->d5s2();

				kern = k[10]; //kD = k[5];
			}
			break;
		case 14:
			{
				_ppb = ppb[14];
				_psb = psb[14];
				_pb = pb[14];

				_v = v[7];

				L6input = p6->d7s1();

				q = a[7];

				dir = 7;

				input = L5A->d7s1();

				kern = k[9];	//kB = k[3];

			}
			break;
		case 15:
			{
				_ppb = ppb[15];
				_psb = psb[15];
				_pb = pb[15];

				_v = v[15];

				L6input = p6->d7s2();

				q = a[15];

				dir = 7;

				input = L5A->d7s2();

				kern = k[11];	//kD = k[7];
			}
			break;
		default:
			throw runtime_error("CLevel5B::compute: invalid direction");
		}

		//Rectify
		rectify(L6input, _ppb, p6->vSize(), p6->hSize());

		#pragma omp barrier

		//Convolve isotropic filter
		if(!(id%2))
			pOs1->conv_sep(_ppb, p6->vSize(), p6->hSize(), _v, _psb); //Convolve for scale 1
		else
			pOs2->conv_sep(_ppb, p6->vSize(), p6->hSize(), _v, _psb); //Convolve for scale 2

		#pragma omp barrier

		//Update all qties
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
				if(!(id%2))
					updateSep(q,dir,input,*kern,pObject,_pb,_psb,_ppb,pObject->_thetas1,L6input);
				else
					updateSep(q,dir,input,*kern,pObject,_pb,_psb,_ppb,pObject->_thetas2,L6input,2);
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
				if(!(id%2))
					update(q,dir,input,*kern,pObject,_pb,_ppb,pObject->_thetas1,L6input);
				else
					update(q,dir,input,*kern,pObject,_pb,_ppb,pObject->_thetas2,L6input,2);
			}
			break;
		default:
			throw runtime_error("CLevel5B::compute: invalid direction");
		}
	}

#endif

#endif

	if(keep_record())
		save();
}

#if NUM_THREADS ==1
double* CLevel5B::update(double* array, unsigned int direction,double* L5AInput, Kernel<double>& k , double theta,double* L6Selection)
{

	static double w0,w1,w2,w3,w4,w5,w6,w7;
	unsigned int offset = _nb_dir - direction;

	w0 = _dg[offset];
	w1 = _dg[offset+1];
	w2 = _dg[offset+2];
	w3 = _dg[offset+3];
	w4 = _dg[offset+4];
	w5 = _dg[offset+5];
	w6 = _dg[offset+6];
	w7 = _dg[offset+7];

	rectify_square(L5AInput,_basePB[0],_vSize,_hSize);
	k.conv_dft(_basePB[0],_baseB[0],_vSize,_hSize);

	double* p = array;
	double* q = _baseB[0];

	double* pT0 = _baseV[0];
	double* pT1 = _baseV[1];
	double* pT2 = _baseV[2];
	double* pT3 = _baseV[3];
	double* pT4 = _baseV[4];
	double* pT5 = _baseV[5];
	double* pT6 = _baseV[6];
	double* pT7 = _baseV[7];

	double* pL6S = L6Selection;

	for(unsigned int i=0;i<_vSize*_hSize;i++,p++,q++,pT0++,pT1++,pT2++,pT3++,pT4++,pT5++,pT6++,pT7++,pL6S++)
	{
		*p = _cst25 * *p +  (1 - *p) * ( ((*q - theta)>0)?(*q - theta):0) * (_cst26 + _cst30*((*pL6S>0)?(*pL6S):0))
			+ (1 + *p) *
			(w0 * *pT0 +
			 w1 * *pT1 +
			 w2 * *pT2 +
			 w3 * *pT3 +
			 w4 * *pT4 +
			 w5 * *pT5 +
			 w6 * *pT6 +
			 w7 * *pT7 );
	}

	return array;
}




double* CLevel5B::updateSep(double* array, unsigned int direction,double* L5AInput, Kernel<double>& k, double theta,double* L6Selection)
{

	static double w0,w1,w2,w3,w4,w5,w6,w7;
	unsigned int offset = _nb_dir - direction;

	w0 = _dg[offset];
	w1 = _dg[offset+1];
	w2 = _dg[offset+2];
	w3 = _dg[offset+3];
	w4 = _dg[offset+4];
	w5 = _dg[offset+5];
	w6 = _dg[offset+6];
	w7 = _dg[offset+7];

	rectify_square(L5AInput,_basePB[0],_vSize,_hSize);
	k.conv_sep(_basePB[0],_vSize,_hSize,_baseB[0],_baseSB[0]);

	double* p = array;
	double* q = _baseB[0];

	double* pT0 = _baseV[0];
	double* pT1 = _baseV[1];
	double* pT2 = _baseV[2];
	double* pT3 = _baseV[3];
	double* pT4 = _baseV[4];
	double* pT5 = _baseV[5];
	double* pT6 = _baseV[6];
	double* pT7 = _baseV[7];

	double* pL6S = L6Selection;

	for(unsigned int i=0;i<_vSize*_hSize;i++,p++,q++,pT0++,pT1++,pT2++,pT3++,pT4++,pT5++,pT6++,pT7++,pL6S++)
	{
		*p = _cst25 * *p + (1 - *p) * ( ((*q - theta)>0)?(*q - theta):0) * (_cst26 + _cst30*((*pL6S>0)?(*pL6S):0))
			+ (1 + *p) *
			(w0 * *pT0 +
			 w1 * *pT1 +
			 w2 * *pT2 +
			 w3 * *pT3 +
			 w4 * *pT4 +
			 w5 * *pT5 +
			 w6 * *pT6 +
			 w7 * *pT7 );
	}

	return array;
}


#else

double* CLevel5B::update(double* array, unsigned int direction, double* L5AInput, Kernel<double>& k,
		CLevel5B* pLayer, double* b, double* preb, double theta, double* L6Selection, int scale)
{
	double w0, w1, w2, w3, w4, w5, w6, w7;
	unsigned int offset = pLayer->_nb_dir - direction;

	w0 = pLayer->_dg[offset];
	w1 = pLayer->_dg[offset+1];
	w2 = pLayer->_dg[offset+2];
	w3 = pLayer->_dg[offset+3];
	w4 = pLayer->_dg[offset+4];
	w5 = pLayer->_dg[offset+5];
	w6 = pLayer->_dg[offset+6];
	w7 = pLayer->_dg[offset+7];

	rectify_square(L5AInput, preb, pLayer->_vSize, pLayer->_hSize);
	k.conv_dft(preb, b, pLayer->_vSize, pLayer->_hSize);

	double* p = array;
	double* q = b;

	double* pT0;
	double* pT1;
	double* pT2;
	double* pT3;
	double* pT4;
	double* pT5;
	double* pT6;
	double* pT7;

#if NUM_THREADS ==16
	if (scale==-1 || scale==1)
	{
#endif
		pT0 = pLayer->_baseV[0];
		pT1 = pLayer->_baseV[1];
		pT2 = pLayer->_baseV[2];
		pT3 = pLayer->_baseV[3];
		pT4 = pLayer->_baseV[4];
		pT5 = pLayer->_baseV[5];
		pT6 = pLayer->_baseV[6];
		pT7 = pLayer->_baseV[7];
#if NUM_THREADS == 16
	}
	else
	{
		pT0 = pLayer->_baseV[8];
		pT1 = pLayer->_baseV[9];
		pT2 = pLayer->_baseV[10];
		pT3 = pLayer->_baseV[11];
		pT4 = pLayer->_baseV[12];
		pT5 = pLayer->_baseV[13];
		pT6 = pLayer->_baseV[14];
		pT7 = pLayer->_baseV[15];
	}
#endif

	double* pL6S = L6Selection;

	for(unsigned int i = 0;i < pLayer->_vSize * pLayer->_hSize; i++, p++, q++, pT0++, pT1++, pT2++, pT3++, pT4++, pT5++, pT6++, pT7++, pL6S++)
	{
		*p = pLayer->_cst25 * *p + (1 - *p) * ( ((*q - theta)>0)?(*q - theta):0) * (pLayer->_cst26 + pLayer->_cst30 * ((*pL6S>0)?(*pL6S):0))
			+ (1 + *p) *
			(w0 * *pT0 +
			 w1 * *pT1 +
			 w2 * *pT2 +
			 w3 * *pT3 +
			 w4 * *pT4 +
			 w5 * *pT5 +
			 w6 * *pT6 +
			 w7 * *pT7 );
	}

	return array;
}

double* CLevel5B::updateSep(double* array, unsigned int direction, double* L5AInput, Kernel<double>& k,
		CLevel5B* pLayer, double* b, double* sb, double* preb,
		double theta, double* L6Selection, int scale)
{

	double w0,w1,w2,w3,w4,w5,w6,w7;
	unsigned int offset = pLayer->_nb_dir - direction;

	w0 = pLayer->_dg[offset];
	w1 = pLayer->_dg[offset+1];
	w2 = pLayer->_dg[offset+2];
	w3 = pLayer->_dg[offset+3];
	w4 = pLayer->_dg[offset+4];
	w5 = pLayer->_dg[offset+5];
	w6 = pLayer->_dg[offset+6];
	w7 = pLayer->_dg[offset+7];

	rectify_square(L5AInput, preb,pLayer->_vSize, pLayer->_hSize);
	k.conv_sep(preb, pLayer->_vSize, pLayer->_hSize, b, sb);

	double* p = array;
	double* q = b;

	double* pT0;
	double* pT1;
	double* pT2;
	double* pT3;
	double* pT4;
	double* pT5;
	double* pT6;
	double* pT7;

#if NUM_THREADS == 16
	if (scale==-1 || scale==1)
	{
#endif

		pT0 = pLayer->_baseV[0];
		pT1 = pLayer->_baseV[1];
		pT2 = pLayer->_baseV[2];
		pT3 = pLayer->_baseV[3];
		pT4 = pLayer->_baseV[4];
		pT5 = pLayer->_baseV[5];
		pT6 = pLayer->_baseV[6];
		pT7 = pLayer->_baseV[7];

#if NUM_THREADS == 16
	}
	else
	{
		pT0 = pLayer->_baseV[8];
		pT1 = pLayer->_baseV[9];
		pT2 = pLayer->_baseV[10];
		pT3 = pLayer->_baseV[11];
		pT4 = pLayer->_baseV[12];
		pT5 = pLayer->_baseV[13];
		pT6 = pLayer->_baseV[14];
		pT7 = pLayer->_baseV[15];
	}
#endif

	double* pL6S = L6Selection;

	for(unsigned int i = 0; i < pLayer->_vSize * pLayer->_hSize; i++, p++, q++, pT0++, pT1++, pT2++, pT3++, pT4++, pT5++, pT6++, pT7++, pL6S++)
	{
		*p = pLayer->_cst25 * *p +  (1 - *p) * ( ((*q - theta)>0)?(*q - theta):0) *(pLayer->_cst26 + pLayer->_cst30 * ((*pL6S>0)?(*pL6S):0))
			+ (1 + *p) *
			(w0 * *pT0 +
			 w1 * *pT1 +
			 w2 * *pT2 +
			 w3 * *pT3 +
			 w4 * *pT4 +
			 w5 * *pT5 +
			 w6 * *pT6 +
			 w7 * *pT7 );
	}
	return array;
}

#endif

void CLevel5B::print_kernels(const string& output_folder, bool recreate, string optName)
{
	//1. Anisotropic kernels
	vector< Kernel<double>* > ks(16);
	ks[0] = &_d0s1K;ks[1] = &_d1s1K;ks[2]  = &_d2s1K;ks[3] = &_d3s1K;ks[4] = &_d0s1K;ks[5] = &_d1s1K;ks[6] = &_d2s1K;ks[7] = &_d3s1K;
	ks[8] = &_d0s2K;ks[9] = &_d1s2K;ks[10]  = &_d2s2K;ks[11] = &_d3s2K;ks[12] = &_d0s2K;ks[13] = &_d1s2K;ks[14] = &_d2s2K;ks[15] = &_d3s2K;

	vector<string> ns(16);
	ns[0] = "D0S1.txt";ns[1] = "D1S1.txt";ns[2] = "D2S1.txt";ns[3] = "D3S1.txt";ns[4] = "D4S1.txt";ns[5] = "D5S1.txt";ns[6] = "D6S1.txt";ns[7] = "D7S1.txt";
	ns[8] = "D0S2.txt";ns[9] = "D1S2.txt";ns[10] = "D2S2.txt";ns[11] = "D3S2.txt";ns[12] = "D4S2.txt";ns[13] = "D5S2.txt";ns[14] = "D6S2.txt";ns[15] = "D7S2.txt";

	vector<string>::iterator n = ns.begin();

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
	string endName = "IsotropicKernelS1";
	base += endName;
	base += optName;
	string endName2 = ".txt";
	base += endName2;

	_p_s1.print_kernel(base, recreate);

	endName = "IsotropicKernelS2";
	base = _level_name;
	base += endName;
	base += optName;
	base += endName2;
	_p_s2.print_kernel(base, recreate);
}
}
