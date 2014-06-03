#include "tinyxml/tinyxml.h"
#include <omp.h>
#include <time.h>
#include <fstream>

#include "visio/KernelLayer.h"
#include "visio/Level6.h"

using namespace std;

namespace visio
{
CLevel6::CLevel6(int iter_per_frame, double dt,
		double A9, double D9, double C9,double B9,
		double A, double sigma,
		double sigmaP,
		unsigned int verticalSize, unsigned int horizontalSize,
		unsigned int totFrames,
		unsigned int nbDirections,
		double V, double sigmaV,
		double VTilde, double sigmaVTilde,
		CLevel5B* pL5,
		Input<int>* pInput, const string& output_folder,
		string name, bool record, Buffer* buffer, bool testKernels)
: KernelLayer(verticalSize, horizontalSize, totFrames, nbDirections,
		sigma, sigma, sigma, sigma, name,
		output_folder, iter_per_frame, dt, record, buffer), _input(pInput)
{

	pL5->set_L6(this);

	//Set constants
	_cst27 = 1-dt*A9;
	_cst28 = -dt*A9*C9;
	_cst29 = dt*A9;
	_cst30 = B9;

	_dg = pL5->dg();

	//Create additional kernel for P
	_P = Kernel<double>(0,sigmaP,sigmaP,_clip_thres);
	_P.weight( (-dt*A9*D9) / (2*PI*sigmaP*sigmaP) );
	_P.set_intrinsic_weight(1 / (2*PI*sigmaP*sigmaP) );

	//SVD for symmetric kernel
	_P.svd();

	//Cross-directional pooling kernel
	vector<double> shape(nbDirections);
	vector<double> shapeTilde(nbDirections);

	_V = vector<double>(2*nbDirections);
	_V_tilde = vector<double>(2*nbDirections);

	for(unsigned int i = 0; i < nbDirections; i++)
	{
		double cDirection = i*2*PI/nbDirections;
		double geodesicDistance = min(cDirection,2*PI - cDirection);
		_V[i] = V * exp(- pow(geodesicDistance/sigmaV,2.0));
		_V_tilde[i] = VTilde * exp(- pow(geodesicDistance/sigmaVTilde,2.0));
	}

	copy(_V.begin(), _V.begin() + nbDirections, _V.begin() + nbDirections);
	copy(_V_tilde.begin(), _V_tilde.begin() + nbDirections, _V_tilde.begin() + nbDirections);

	setBaseSB();
	setBasePB();
	setBaseV();
	setBaseRL5();
	setBaseRL6();

	if(testKernels)
	{
		time_kernel(verticalSize, horizontalSize);
		time_kernel(verticalSize, horizontalSize, &_P, "P");
	}
}

void CLevel6::setBaseSB()
{
#if NUM_THREADS == 1
	_baseSB.push_back(_buf->_buffer0);
#elif NUM_THREADS == 4
	_baseSB.push_back(_buf->_buffer0);
	_baseSB.push_back(_buf->_buffer1);
	_baseSB.push_back(_buf->_buffer2);
	_baseSB.push_back(_buf->_buffer3);
#elif NUM_THREADS == 16
	_baseSB.push_back(_buf->_buffer0);
	_baseSB.push_back(_buf->_buffer1);
	_baseSB.push_back(_buf->_buffer2);
	_baseSB.push_back(_buf->_buffer3);
	_baseSB.push_back(_buf->_buffer4);
	_baseSB.push_back(_buf->_buffer5);
	_baseSB.push_back(_buf->_buffer6);
	_baseSB.push_back(_buf->_buffer7);
	_baseSB.push_back(_buf->_buffer8);
	_baseSB.push_back(_buf->_buffer9);
	_baseSB.push_back(_buf->_buffer10);
	_baseSB.push_back(_buf->_buffer11);
	_baseSB.push_back(_buf->_buffer12);
	_baseSB.push_back(_buf->_buffer13);
	_baseSB.push_back(_buf->_buffer14);
	_baseSB.push_back(_buf->_buffer15);
#endif
}

void CLevel6::setBasePB()
{
#if NUM_THREADS == 4
	_basePB.push_back(_buf->_buffer12);
	_basePB.push_back(_buf->_buffer13);
	_basePB.push_back(_buf->_buffer14);
	_basePB.push_back(_buf->_buffer15);
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

void CLevel6::setBaseV()
{
#if NUM_THREADS == 1
	_baseV.push_back(_buf->_buffer1);
	_baseV.push_back(_buf->_buffer2);
	_baseV.push_back(_buf->_buffer3);
	_baseV.push_back(_buf->_buffer4);
	_baseV.push_back(_buf->_buffer5);
	_baseV.push_back(_buf->_buffer6);
	_baseV.push_back(_buf->_buffer7);
	_baseV.push_back(_buf->_buffer8);
#elif NUM_THREADS == 4
	_baseV.push_back(_buf->_buffer4);
	_baseV.push_back(_buf->_buffer5);
	_baseV.push_back(_buf->_buffer6);
	_baseV.push_back(_buf->_buffer7);
	_baseV.push_back(_buf->_buffer8);
	_baseV.push_back(_buf->_buffer9);
	_baseV.push_back(_buf->_buffer10);
	_baseV.push_back(_buf->_buffer11);
#elif NUM_THREADS == 16
	_baseV.push_back(_buf->_buffer16);
	_baseV.push_back(_buf->_buffer17);
	_baseV.push_back(_buf->_buffer18);
	_baseV.push_back(_buf->_buffer19);
	_baseV.push_back(_buf->_buffer20);
	_baseV.push_back(_buf->_buffer21);
	_baseV.push_back(_buf->_buffer22);
	_baseV.push_back(_buf->_buffer23);
	_baseV.push_back(_buf->_buffer24);
	_baseV.push_back(_buf->_buffer25);
	_baseV.push_back(_buf->_buffer26);
	_baseV.push_back(_buf->_buffer27);
	_baseV.push_back(_buf->_buffer28);
	_baseV.push_back(_buf->_buffer29);
	_baseV.push_back(_buf->_buffer30);
	_baseV.push_back(_buf->_buffer31);
#endif
}

void CLevel6::setBaseRL5()
{
#if NUM_THREADS == 1
	_baseRL5.push_back(_buf->_buffer9);
	_baseRL5.push_back(_buf->_buffer10);
	_baseRL5.push_back(_buf->_buffer11);
	_baseRL5.push_back(_buf->_buffer12);
	_baseRL5.push_back(_buf->_buffer13);
	_baseRL5.push_back(_buf->_buffer14);
	_baseRL5.push_back(_buf->_buffer15);
	_baseRL5.push_back(_buf->_buffer16);
#elif NUM_THREADS == 4
	_baseRL5.push_back(_buf->_buffer16);
	_baseRL5.push_back(_buf->_buffer17);
	_baseRL5.push_back(_buf->_buffer18);
	_baseRL5.push_back(_buf->_buffer19);
	_baseRL5.push_back(_buf->_buffer20);
	_baseRL5.push_back(_buf->_buffer21);
	_baseRL5.push_back(_buf->_buffer22);
	_baseRL5.push_back(_buf->_buffer23);
#elif NUM_THREADS == 16
	_baseRL5.push_back(_buf->_buffer48);
	_baseRL5.push_back(_buf->_buffer49);
	_baseRL5.push_back(_buf->_buffer50);
	_baseRL5.push_back(_buf->_buffer51);
	_baseRL5.push_back(_buf->_buffer52);
	_baseRL5.push_back(_buf->_buffer53);
	_baseRL5.push_back(_buf->_buffer54);
	_baseRL5.push_back(_buf->_buffer55);
	_baseRL5.push_back(_buf->_buffer56);
	_baseRL5.push_back(_buf->_buffer57);
	_baseRL5.push_back(_buf->_buffer58);
	_baseRL5.push_back(_buf->_buffer59);
	_baseRL5.push_back(_buf->_buffer60);
	_baseRL5.push_back(_buf->_buffer61);
	_baseRL5.push_back(_buf->_buffer62);
	_baseRL5.push_back(_buf->_buffer63);
#endif
}

void CLevel6::setBaseRL6()
{
#if NUM_THREADS == 1
	_baseRL6.push_back(_buf->_buffer17);
	_baseRL6.push_back(_buf->_buffer18);
	_baseRL6.push_back(_buf->_buffer19);
	_baseRL6.push_back(_buf->_buffer20);
	_baseRL6.push_back(_buf->_buffer21);
	_baseRL6.push_back(_buf->_buffer22);
	_baseRL6.push_back(_buf->_buffer23);
	_baseRL6.push_back(_buf->_buffer24);
#elif NUM_THREADS == 4
	_baseRL6.push_back(_buf->_buffer24);
	_baseRL6.push_back(_buf->_buffer25);
	_baseRL6.push_back(_buf->_buffer26);
	_baseRL6.push_back(_buf->_buffer27);
	_baseRL6.push_back(_buf->_buffer28);
	_baseRL6.push_back(_buf->_buffer29);
	_baseRL6.push_back(_buf->_buffer30);
	_baseRL6.push_back(_buf->_buffer31);
#elif NUM_THREADS == 16
	_baseRL6.push_back(_buf->_buffer64);
	_baseRL6.push_back(_buf->_buffer65);
	_baseRL6.push_back(_buf->_buffer66);
	_baseRL6.push_back(_buf->_buffer67);
	_baseRL6.push_back(_buf->_buffer68);
	_baseRL6.push_back(_buf->_buffer69);
	_baseRL6.push_back(_buf->_buffer70);
	_baseRL6.push_back(_buf->_buffer71);
#endif
}

void CLevel6::compute(void* object)
{
	CLevel5B* L5B = static_cast<CLevel5B*>(object);

#if NUM_THREADS==1

	//1. Convolve T once for scale 1
	_P.conv_sepRectified(m_d0s1,_vSize,_hSize,_baseV[0],_baseSB[0]);
	_P.conv_sepRectified(m_d1s1,_vSize,_hSize,_baseV[1],_baseSB[0]);
	_P.conv_sepRectified(m_d2s1,_vSize,_hSize,_baseV[2],_baseSB[0]);
	_P.conv_sepRectified(m_d3s1,_vSize,_hSize,_baseV[3],_baseSB[0]);
	_P.conv_sepRectified(m_d4s1,_vSize,_hSize,_baseV[4],_baseSB[0]);
	_P.conv_sepRectified(m_d5s1,_vSize,_hSize,_baseV[5],_baseSB[0]);
	_P.conv_sepRectified(m_d6s1,_vSize,_hSize,_baseV[6],_baseSB[0]);
	_P.conv_sepRectified(m_d7s1,_vSize,_hSize,_baseV[7],_baseSB[0]);

	rectify(L5B->array("D0S1"),_baseRL5[0],_vSize,_hSize);
	rectify(L5B->array("D1S1"),_baseRL5[1],_vSize,_hSize);
	rectify(L5B->array("D2S1"),_baseRL5[2],_vSize,_hSize);
	rectify(L5B->array("D3S1"),_baseRL5[3],_vSize,_hSize);
	rectify(L5B->array("D4S1"),_baseRL5[4],_vSize,_hSize);
	rectify(L5B->array("D5S1"),_baseRL5[5],_vSize,_hSize);
	rectify(L5B->array("D6S1"),_baseRL5[6],_vSize,_hSize);
	rectify(L5B->array("D7S1"),_baseRL5[7],_vSize,_hSize);

	//Update at scale 1 (no depth inhibition)
	if(_input->attend_direction(0))
		update_full_no_depth(m_d0s1,0,this);
	else
		update_noO_noDepth_gaussian(m_d0s1,0,this);

	if(_input->attend_direction(1))
		update_full_no_depth(m_d1s1,1,this);
	else
		update_noO_noDepth_gaussian(m_d1s1,1,this);

	if(_input->attend_direction(2))
		update_full_no_depth(m_d2s1,2,this);
	else
		update_noO_noDepth_gaussian(m_d2s1,2,this);

	if(_input->attend_direction(3))
		update_full_no_depth(m_d3s1,3,this);
	else
		update_noO_noDepth_gaussian(m_d3s1,3,this);

	if(_input->attend_direction(4))
		update_full_no_depth(m_d4s1,4,this);
	else
		update_noO_noDepth_gaussian(m_d4s1,4,this);

	if(_input->attend_direction(5))
		update_full_no_depth(m_d5s1,5,this);
	else
		update_noO_noDepth_gaussian(m_d5s1,5,this);

	if(_input->attend_direction(6))
		update_full_no_depth(m_d6s1,6,this);
	else
		update_noO_noDepth_gaussian(m_d6s1,6,this);

	if(_input->attend_direction(7))
		update_full_no_depth(m_d7s1,7,this);
	else
		update_noO_noDepth_gaussian(m_d7s1,7,this);

	//1. Convolve T once for scale 2
	_P.conv_sepRectified(m_d0s2,_vSize,_hSize,_baseV[0],_baseSB[0]);
	_P.conv_sepRectified(m_d1s2,_vSize,_hSize,_baseV[1],_baseSB[0]);
	_P.conv_sepRectified(m_d2s2,_vSize,_hSize,_baseV[2],_baseSB[0]);
	_P.conv_sepRectified(m_d3s2,_vSize,_hSize,_baseV[3],_baseSB[0]);
	_P.conv_sepRectified(m_d4s2,_vSize,_hSize,_baseV[4],_baseSB[0]);
	_P.conv_sepRectified(m_d5s2,_vSize,_hSize,_baseV[5],_baseSB[0]);
	_P.conv_sepRectified(m_d6s2,_vSize,_hSize,_baseV[6],_baseSB[0]);
	_P.conv_sepRectified(m_d7s2,_vSize,_hSize,_baseV[7],_baseSB[0]);

	//Rectify the activity at scale 1, for future use in depth suppression
	rectify(m_d0s1,_baseRL6[0],_vSize,_hSize);
	rectify(m_d1s1,_baseRL6[1],_vSize,_hSize);
	rectify(m_d2s1,_baseRL6[2],_vSize,_hSize);
	rectify(m_d3s1,_baseRL6[3],_vSize,_hSize);
	rectify(m_d4s1,_baseRL6[4],_vSize,_hSize);
	rectify(m_d5s1,_baseRL6[5],_vSize,_hSize);
	rectify(m_d6s1,_baseRL6[6],_vSize,_hSize);
	rectify(m_d7s1,_baseRL6[7],_vSize,_hSize);

	//Update at scale 2	(depth inhibition from plane 1)
	rectify(L5B->array("D0S2"),_baseRL5[0],_vSize,_hSize);
	rectify(L5B->array("D1S2"),_baseRL5[1],_vSize,_hSize);
	rectify(L5B->array("D2S2"),_baseRL5[2],_vSize,_hSize);
	rectify(L5B->array("D3S2"),_baseRL5[3],_vSize,_hSize);
	rectify(L5B->array("D4S2"),_baseRL5[4],_vSize,_hSize);
	rectify(L5B->array("D5S2"),_baseRL5[5],_vSize,_hSize);
	rectify(L5B->array("D6S2"),_baseRL5[6],_vSize,_hSize);
	rectify(L5B->array("D7S2"),_baseRL5[7],_vSize,_hSize);

	update_noO_depth_gaussian(m_d0s2,0,this);
	update_noO_depth_gaussian(m_d1s2,1,this);
	update_noO_depth_gaussian(m_d2s2,2,this);
	update_noO_depth_gaussian(m_d3s2,3,this);
	update_noO_depth_gaussian(m_d4s2,4,this);
	update_noO_depth_gaussian(m_d5s2,5,this);
	update_noO_depth_gaussian(m_d6s2,6,this);
	update_noO_depth_gaussian(m_d7s2,7,this);

#else

	Kernel<double>* pP = &_P;

	vector<double*> pb   = _baseB;
	vector<double*> psb  = _baseSB;
	vector<double*> ppb  = _basePB;
	vector<double*> v    = _baseV;
	vector<double*> prL5 = _baseRL5;
	vector<double*> prL6 = _baseRL6;

	auto k = _k;
	vector<double*> a = _a;

	CLevel6* pObject = this;
	Input<int>* pInputLayer = _input;

#if NUM_THREADS == 4

	omp_set_num_threads(4);

	//Parrallelize cvd computation
	#pragma omp parallel shared(pP,a,v,psb,ppb,pObject,L5B,prL5,prL6,pInputLayer)
	{
		int id = omp_get_thread_num();

		double* _ppb;
		double* _psb;
		double* _vA;
		double* _vB;
		double* input1;
		double* input2;
		double* input3;
		double* input4;
		double* q1;
		double* q2;
		double* q3;
		double* q4;

		Kernel<double>* kA;
		Kernel<double>* kB;

		unsigned int dirA;
		unsigned int dirB;

		double* hA;
		double* hB;

		double* r5A;
		double* r5B;
		double* r6A;
		double* r6B;

		switch(id)
		{
		case 0:
			{
				q1   = pObject->d0s1();
				q2   = pObject->d1s1();
				q3   = pObject->d0s2();
				q4   = pObject->d1s2();

				input1 = L5B->d0s1();
				input2 = L5B->d1s1();
				input3 = L5B->d0s2();
				input4 = L5B->d1s2();

				_ppb = ppb[0];
				_psb = psb[0];
				_vA  = v[0];
				_vB  = v[1];

				kA = k[0];
				kB = k[1];

				dirA = 0;
				dirB = 1;

				hA = pObject->d0s1();
				hB = pObject->d1s1();

				r5A = prL5[0];
				r5B = prL5[1];
				r6A = prL6[0];
				r6B = prL6[1];
			}
			break;
		case 1:
			{
				q1   = pObject->d2s1();
				q2   = pObject->d3s1();
				q3   = pObject->d2s2();
				q4   = pObject->d3s2();

				input1 = L5B->d2s1();
				input2 = L5B->d3s1();
				input3 = L5B->d2s2();
				input4 = L5B->d3s2();

				_ppb = ppb[1];
				_psb = psb[1];
				_vA  = v[2];
				_vB  = v[3];

				kA = k[2];
				kB = k[3];

				dirA = 2;
				dirB = 3;

				hA = pObject->d2s1();
				hB = pObject->d3s1();

				r5A = prL5[2];
				r5B = prL5[3];
				r6A = prL6[2];
				r6B = prL6[3];
			}
			break;
		case 2:
			{
				q1   = pObject->d4s1();
				q2   = pObject->d5s1();
				q3   = pObject->d4s2();
				q4   = pObject->d5s2();

				input1 = L5B->d4s1();
				input2 = L5B->d5s1();
				input3 = L5B->d4s2();
				input4 = L5B->d5s2();

				_ppb = ppb[2];
				_psb = psb[2];
				_vA  = v[4];
				_vB  = v[5];

				kA = k[0];
				kB = k[8];	//kB = k[1];

				dirA = 4;
				dirB = 5;

				hA = pObject->d4s1();
				hB = pObject->d5s1();

				r5A = prL5[4];
				r5B = prL5[5];
				r6A = prL6[4];
				r6B = prL6[5];
			}
			break;
		case 3:
			{
				q1   = pObject->d6s1();
				q2   = pObject->d7s1();
				q3   = pObject->d6s2();
				q4   = pObject->d7s2();

				input1 = L5B->d6s1();
				input2 = L5B->d7s1();
				input3 = L5B->d6s2();
				input4 = L5B->d7s2();

				_ppb = ppb[3];
				_psb = psb[3];
				_vA  = v[6];
				_vB  = v[7];

				kA = k[2];
				kB = k[9];	//kB = k[3];

				dirA = 6;
				dirB = 7;

				hA = pObject->d6s1();
				hB = pObject->d7s1();

				r5A = prL5[6];
				r5B = prL5[7];
				r6A = prL6[6];
				r6B = prL6[7];
			}
		}

		rectify(q1,_ppb,pObject->vSize(),pObject->hSize());
		pP->conv_sep(_ppb,pObject->vSize(),pObject->hSize(),_vA,_psb);

		rectify(q2,_ppb,pObject->vSize(),pObject->hSize());
		pP->conv_sep(_ppb,pObject->vSize(),pObject->hSize(),_vB,_psb);

		//Rectify Layer 5B input
		rectify(input1,r5A,pObject->vSize(),pObject->hSize());
		rectify(input2,r5B,pObject->vSize(),pObject->hSize());

		#pragma omp barrier

		if(pInputLayer->attend_direction(dirA))		update_full_no_depth(q1,dirA,pObject);
		else									update_noO_noDepth_gaussian(q1,dirA,pObject);

		if(pInputLayer->attend_direction(dirB))		update_full_no_depth(q2,dirB,pObject);
		else									update_noO_noDepth_gaussian(q2,dirB,pObject);

		#pragma omp barrier

		rectify(q3,_ppb,pObject->vSize(),pObject->hSize());
		pP->conv_sep(_ppb,pObject->vSize(),pObject->hSize(),_vA,_psb);

		rectify(q4,_ppb,pObject->vSize(),pObject->hSize());
		pP->conv_sep(_ppb,pObject->vSize(),pObject->hSize(),_vB,_psb);

		//Rectify Layer 5B input
		rectify(input3,r5A,pObject->vSize(),pObject->hSize());
		rectify(input4,r5B,pObject->vSize(),pObject->hSize());

		//Rectify Layer 6, scale 1, signals
		rectify(q1,r6A,pObject->vSize(),pObject->hSize());
		rectify(q2,r6B,pObject->vSize(),pObject->hSize());

		#pragma omp barrier

		update_noO_depth_gaussian(q3,dirA,pObject);
		update_noO_depth_gaussian(q4,dirB,pObject);

		}

#elif NUM_THREADS == 16

	omp_set_num_threads(16);

	//Parrallelize cvd computation
	#pragma omp parallel shared(pP,a,v,psb,ppb,pObject,L5B,prL5,prL6,pInputLayer)
	{
		double* _ppb = nullptr;
		double* _psb = nullptr;
		double* _v = nullptr;

		double* input = nullptr;
		double* q1 = nullptr;
		double* q2 = nullptr;

		unsigned int dir;

		double* r5 = nullptr;
		double* r6 = nullptr;

		int id = omp_get_thread_num();
		switch(id)
		{
		case 0:
			{
				q1   = pObject->d0s1();
			    q2   = pObject->d0s2();

				input = L5B->d0s1();

				_ppb = ppb[0];
				_psb = psb[0];
				_v  = v[0];

				dir = 0;

				r5 = prL5[0];
				r6 = prL6[0];
			}
			break;
		case 1:
			{
				q1   = pObject->d0s2();

				input = L5B->d0s2();

				_ppb = ppb[1];
				_psb = psb[1];
				_v   = v[8];

				dir = 0;

				r5 = prL5[8];
			}
			break;
		case 2:
			{
				q1   = pObject->d2s1();
				q2   = pObject->d2s2();

				input = L5B->d2s1();

				_ppb = ppb[2];
				_psb = psb[2];
				_v   = v[2];

				dir = 2;

				r5 = prL5[2];
				r6 = prL6[2];
			}
			break;
		case 3:
			{
				q1   = pObject->d2s2();

				input = L5B->d2s2();

				_ppb = ppb[3];
				_psb = psb[3];
				_v   = v[10];

				dir = 2;

				r5 = prL5[10];
			}
			break;
		case 4:
			{
				q1   = pObject->d4s1();
				q2   = pObject->d4s2();

				input = L5B->d4s1();

				_ppb = ppb[4];
				_psb = psb[4];
				_v   = v[4];

				dir = 4;

				r5 = prL5[4];
				r6 = prL6[4];
			}
			break;
		case 5:
			{
				q1   = pObject->d4s2();

				input = L5B->d4s2();

				_ppb = ppb[5];
				_psb = psb[5];
				_v   = v[12];

				dir = 4;

				r5 = prL5[12];
			}
			break;
		case 6:
			{
				q1   = pObject->d6s1();
				q2   = pObject->d6s2();

				input = L5B->d6s1();

				_ppb = ppb[6];
				_psb = psb[6];
				_v   = v[6];

				dir = 6;

				r5 = prL5[6];
				r6 = prL6[6];
			}
			break;
		case 7:
			{
				q1   = pObject->d6s2();

				input = L5B->d6s2();

				_ppb = ppb[7];
				_psb = psb[7];
				_v   = v[14];

				dir = 6;

				r5 = prL5[14];
			}
			break;
		case 8:
			{
				q1   = pObject->d1s1();
				q2   = pObject->d1s2();

				input = L5B->d1s1();

				_ppb = ppb[8];
				_psb = psb[8];

				_v   = v[1];

				dir = 1;

				r5 = prL5[1];
				r6 = prL6[1];
			}
			break;
		case 9:
			{
				q1   = pObject->d1s2();

				input = L5B->d1s2();

				_ppb = ppb[9];
				_psb = psb[9];
				_v   = v[9];

				dir = 1;

				r5 = prL5[9];
			}
			break;
		case 10:
			{
				q1   = pObject->d3s1();
				q2   = pObject->d3s2();

				input = L5B->d3s1();

				_ppb = ppb[10];
				_psb = psb[10];
				_v   = v[3];

				dir = 3;

				r5 = prL5[3];
				r6 = prL6[3];
			}
			break;
		case 11:
			{
				q1   = pObject->d3s2();

				input = L5B->d3s2();

				_ppb = ppb[11];
				_psb = psb[11];
				_v   = v[11];

				dir = 3;

				r5 = prL5[11];
			}
			break;
		case 12:
			{
				q1   = pObject->d5s1();
				q2   = pObject->d5s2();

				input = L5B->d5s1();

				_ppb = ppb[12];
				_psb = psb[12];
				_v   = v[5];

				dir = 5;

				r5 = prL5[5];
				r6 = prL6[5];
			}
			break;
		case 13:
			{
				q1   = pObject->d5s2();

				input = L5B->d5s2();

				_ppb = ppb[13];
				_psb = psb[13];
				_v   = v[13];

				dir = 5;

				r5 = prL5[13];
			}
			break;
		case 14:
			{
				q1   = pObject->d7s1();
				q2   = pObject->d7s2();

				input = L5B->d7s1();

				_ppb = ppb[14];
				_psb = psb[14];
				_v   = v[7];

				dir = 7;

				r5 = prL5[7];
				r6 = prL6[7];
			}
			break;
		case 15:
			{
				q1   = pObject->d7s2();

				input = L5B->d7s2();

				_ppb = ppb[15];
				_psb = psb[15];

				_v   = v[15];

				dir = 7;

				r5 = prL5[15];
			}
			break;
		default:
			throw runtime_error("CLevel6::compute: invalid direction");
		}

		//Rectify Level 6 activity
		rectify(q1, _ppb, pObject->vSize(), pObject->hSize());

		//Convolve with isotropic kernel
		pP->conv_sep(_ppb, pObject->vSize(), pObject->hSize(), _v, _psb);

		//Rectify Level 5B input
		rectify(input, r5, pObject->vSize(), pObject->hSize());

		#pragma omp barrier

		if(!(id%2))
		{
			//If is an even numbered thread, first update scale 1 (large, near) of level 6
			if(pInputLayer->attend_direction(dir))
				update_full_no_depth(q1, dir, pObject);
			else
				update_noO_noDepth_gaussian(q1, dir, pObject);

			#pragma omp barrier

			//Rectify Layer 6, scale 1, signals
			rectify(q1, r6, pObject->vSize(), pObject->hSize());

			#pragma omp barrier

			//Update scale 2 (far, small) of level 6, with depth suppression
			update_noO_depth_gaussian(q2, dir, pObject);
		}
		else
		{
			//Otherwise, do nothing
			#pragma omp barrier

			#pragma omp barrier
		}
	}

#endif

#endif

	if(keep_record())
		save();
}

double* CLevel6::update_full_no_depth(double* array, unsigned int direction, CLevel6* pLayer)
{
	double* p = array;

	double w0,w1,w2,w3,w4,w5,w6,w7;
	double v0,v1,v2,v3,v4,v5,v6,v7;

	unsigned int offset = pLayer->_nb_dir - direction;

	w0 = pLayer->_dg[offset];
	w1 = pLayer->_dg[offset+1];
	w2 = pLayer->_dg[offset+2];
	w3 = pLayer->_dg[offset+3];
	w4 = pLayer->_dg[offset+4];
	w5 = pLayer->_dg[offset+5];
	w6 = pLayer->_dg[offset+6];
	w7 = pLayer->_dg[offset+7];

	v0 = pLayer->_V[offset];
	v1 = pLayer->_V[offset+1];
	v2 = pLayer->_V[offset+2];
	v3 = pLayer->_V[offset+3];
	v4 = pLayer->_V[offset+4];
	v5 = pLayer->_V[offset+5];
	v6 = pLayer->_V[offset+6];
	v7 = pLayer->_V[offset+7];

	double* pT0 = pLayer->_baseV[0];
	double* pT1 = pLayer->_baseV[1];
	double* pT2 = pLayer->_baseV[2];
	double* pT3 = pLayer->_baseV[3];
	double* pT4 = pLayer->_baseV[4];
	double* pT5 = pLayer->_baseV[5];
	double* pT6 = pLayer->_baseV[6];
	double* pT7 = pLayer->_baseV[7];

	double* pL5d0 = pLayer->_baseRL5[0];
	double* pL5d1 = pLayer->_baseRL5[1];
	double* pL5d2 = pLayer->_baseRL5[2];
	double* pL5d3 = pLayer->_baseRL5[3];
	double* pL5d4 = pLayer->_baseRL5[4];
	double* pL5d5 = pLayer->_baseRL5[5];
	double* pL5d6 = pLayer->_baseRL5[6];
	double* pL5d7 = pLayer->_baseRL5[7];

	//Attentional parameters
	vector<double> attParam;
	_input->attention(direction, attParam);
	unsigned int x, y;

	for(unsigned int i = 0;i < pLayer->_vSize * pLayer->_hSize; i++,p++,pT0++,pT1++,pT2++,pT3++,pT4++,pT5++,pT6++,pT7++,
														   pL5d0++,pL5d1++,pL5d2++,pL5d3++,pL5d4++,pL5d5++,pL5d6++,pL5d7++)
	{
		y = i/pLayer->_hSize + 1;  //Get current coordinates
		x = i%pLayer->_hSize + 1;	//Add (1,1) since reference frame for attentional kernel starts in (1,1), whereas it starts in (0,0) for index i

		*p = pLayer->_cst27* *p + (1 - *p) * pLayer->_cst29 *
			(
			 v0 * *pL5d0 +
			 v1 * *pL5d1 +
			 v2 * *pL5d2 +
			 v3 * *pL5d3 +
			 v4 * *pL5d4 +
			 v5 * *pL5d5 +
			 v6 * *pL5d6 +
			 v7 * *pL5d7
			 ) *
			 (1 + top_down_weight(x, y, attParam[0], attParam[1], attParam[2], attParam[3]))
			 +
			 (_cst30 + *p)	//No need to care about - sign,  dt or A9, since already weighted in the m_P kernel earlier
			 *
			 (
			 w0 * *pT0 +
			 w1 * *pT1 +
			 w2 * *pT2 +
			 w3 * *pT3 +
			 w4 * *pT4 +
			 w5 * *pT5 +
			 w6 * *pT6 +
			 w7 * *pT7);
	}
	return array;
}


double* CLevel6::update_noO_depth(double* array, unsigned int direction, double* L5BInput,
		double* op, CLevel6* pLayer)
{
	double* p = array;
	double* q = L5BInput;
	double* r = op;

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

#if NUM_THREADS == 1 || NUM_THREADS == 4
	double* pT0 = pLayer->_baseV[0];
	double* pT1 = pLayer->_baseV[1];
	double* pT2 = pLayer->_baseV[2];
	double* pT3 = pLayer->_baseV[3];
	double* pT4 = pLayer->_baseV[4];
	double* pT5 = pLayer->_baseV[5];
	double* pT6 = pLayer->_baseV[6];
	double* pT7 = pLayer->_baseV[7];
#elif NUM_THREADS == 16
	double* pT0 = pLayer->_baseV[8];
	double* pT1 = pLayer->_baseV[9];
	double* pT2 = pLayer->_baseV[10];
	double* pT3 = pLayer->_baseV[11];
	double* pT4 = pLayer->_baseV[12];
	double* pT5 = pLayer->_baseV[13];
	double* pT6 = pLayer->_baseV[14];
	double* pT7 = pLayer->_baseV[15];
#endif

	for(unsigned int i = 0;i < pLayer->_vSize * pLayer->_hSize; i++,p++,q++,pT0++,pT1++,pT2++,pT3++,pT4++,pT5++,pT6++,pT7++,r++)
	{
		*p = pLayer->_cst27* *p + (1 - *p) * ( (*q>0)?*q:0 ) * pLayer->_cst29 +
			 (_cst30 + *p)	//No need to care about - sign,  dt or A9, since already weighted in the m_P kernel earlier
			 *
			 (w0 * *pT0 +
			 w1 * *pT1 +
			 w2 * *pT2 +
			 w3 * *pT3 +
			 w4 * *pT4 +
			 w5 * *pT5 +
			 w6 * *pT6 +
			 w7 * *pT7
			 + pLayer->_cst28 * *r);
	}
	return array;
}
double* CLevel6::update_noO_noDepth(double* array, unsigned int direction, double* L5BInput, CLevel6* pLayer)
{
	double* p = array;
	double* q = L5BInput;

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

	double* pT0 = pLayer->_baseV[0];
	double* pT1 = pLayer->_baseV[1];
	double* pT2 = pLayer->_baseV[2];
	double* pT3 = pLayer->_baseV[3];
	double* pT4 = pLayer->_baseV[4];
	double* pT5 = pLayer->_baseV[5];
	double* pT6 = pLayer->_baseV[6];
	double* pT7 = pLayer->_baseV[7];

	for(unsigned int i = 0; i < pLayer->_vSize * pLayer->_hSize; i++,p++,q++,pT0++,pT1++,pT2++,pT3++,pT4++,pT5++,pT6++,pT7++)
	{
		*p = pLayer->_cst27* *p + (1 - *p) * ( (*q>0)?*q:0 ) * pLayer->_cst29 +
			 (_cst30 + *p)	//No need to care about - sign,  dt or A9, since already weighted in the m_P kernel earlier
			 *
			 (w0 * *pT0 +
			 w1 * *pT1 +
			 w2 * *pT2 +
			 w3 * *pT3 +
			 w4 * *pT4 +
			 w5 * *pT5 +
			 w6 * *pT6 +
			 w7 * *pT7);
	}


	return array;
}

double* CLevel6::update_noO_depth_gaussian(double* array, unsigned int direction, CLevel6* pLayer)
{
	double* p = array;

	double w0,w1,w2,w3,w4,w5,w6,w7;
	double v0,v1,v2,v3,v4,v5,v6,v7;
	double vT0,vT1,vT2,vT3,vT4,vT5,vT6,vT7;

	unsigned int offset = pLayer->_nb_dir - direction;

	w0 = pLayer->_dg[offset];
	w1 = pLayer->_dg[offset+1];
	w2 = pLayer->_dg[offset+2];
	w3 = pLayer->_dg[offset+3];
	w4 = pLayer->_dg[offset+4];
	w5 = pLayer->_dg[offset+5];
	w6 = pLayer->_dg[offset+6];
	w7 = pLayer->_dg[offset+7];

	v0 = pLayer->_V[offset];
	v1 = pLayer->_V[offset+1];
	v2 = pLayer->_V[offset+2];
	v3 = pLayer->_V[offset+3];
	v4 = pLayer->_V[offset+4];
	v5 = pLayer->_V[offset+5];
	v6 = pLayer->_V[offset+6];
	v7 = pLayer->_V[offset+7];

	vT0 = pLayer->_V_tilde[offset];
	vT1 = pLayer->_V_tilde[offset+1];
	vT2 = pLayer->_V_tilde[offset+2];
	vT3 = pLayer->_V_tilde[offset+3];
	vT4 = pLayer->_V_tilde[offset+4];
	vT5 = pLayer->_V_tilde[offset+5];
	vT6 = pLayer->_V_tilde[offset+6];
	vT7 = pLayer->_V_tilde[offset+7];

#if NUM_THREADS == 1 || NUM_THREADS == 4

	double* pT0 = pLayer->_baseV[0];
	double* pT1 = pLayer->_baseV[1];
	double* pT2 = pLayer->_baseV[2];
	double* pT3 = pLayer->_baseV[3];
	double* pT4 = pLayer->_baseV[4];
	double* pT5 = pLayer->_baseV[5];
	double* pT6 = pLayer->_baseV[6];
	double* pT7 = pLayer->_baseV[7];

	double* pL5d0 = pLayer->_baseRL5[0];
	double* pL5d1 = pLayer->_baseRL5[1];
	double* pL5d2 = pLayer->_baseRL5[2];
	double* pL5d3 = pLayer->_baseRL5[3];
	double* pL5d4 = pLayer->_baseRL5[4];
	double* pL5d5 = pLayer->_baseRL5[5];
	double* pL5d6 = pLayer->_baseRL5[6];
	double* pL5d7 = pLayer->_baseRL5[7];

#elif NUM_THREADS == 16

	double* pT0 = pLayer->_baseV[8];
	double* pT1 = pLayer->_baseV[9];
	double* pT2 = pLayer->_baseV[10];
	double* pT3 = pLayer->_baseV[11];
	double* pT4 = pLayer->_baseV[12];
	double* pT5 = pLayer->_baseV[13];
	double* pT6 = pLayer->_baseV[14];
	double* pT7 = pLayer->_baseV[15];

	double* pL5d0 = pLayer->_baseRL5[8];
	double* pL5d1 = pLayer->_baseRL5[9];
	double* pL5d2 = pLayer->_baseRL5[10];
	double* pL5d3 = pLayer->_baseRL5[11];
	double* pL5d4 = pLayer->_baseRL5[12];
	double* pL5d5 = pLayer->_baseRL5[13];
	double* pL5d6 = pLayer->_baseRL5[14];
	double* pL5d7 = pLayer->_baseRL5[15];

#endif

	double* pL6d0s1 = pLayer->_baseRL6[0];
	double* pL6d1s1 = pLayer->_baseRL6[1];
	double* pL6d2s1 = pLayer->_baseRL6[2];
	double* pL6d3s1 = pLayer->_baseRL6[3];
	double* pL6d4s1 = pLayer->_baseRL6[4];
	double* pL6d5s1 = pLayer->_baseRL6[5];
	double* pL6d6s1 = pLayer->_baseRL6[6];
	double* pL6d7s1 = pLayer->_baseRL6[7];

	for(unsigned int i = 0; i < pLayer->_vSize * pLayer->_hSize; i++,p++,pT0++,pT1++,pT2++,pT3++,pT4++,pT5++,pT6++,pT7++,
														   pL5d0++,pL5d1++,pL5d2++,pL5d3++,pL5d4++,pL5d5++,pL5d6++,pL5d7++,
														   pL6d0s1++,pL6d1s1++,pL6d2s1++,pL6d3s1++,pL6d4s1++,pL6d5s1++,pL6d6s1++,pL6d7s1++)
	{
		*p = pLayer->_cst27* *p + (1 - *p) * pLayer->_cst29 *
			(
			 v0 * *pL5d0 +
			 v1 * *pL5d1 +
			 v2 * *pL5d2 +
			 v3 * *pL5d3 +
			 v4 * *pL5d4 +
			 v5 * *pL5d5 +
			 v6 * *pL5d6 +
			 v7 * *pL5d7
			 ) +
			 (_cst30 + *p)//No need to care about - sign,  dt or A9, since already weighted in the m_P kernel earlier, and -dtA9C9 is in m_cst28
			 *
			 (w0 * *pT0 +
			 w1 * *pT1 +
			 w2 * *pT2 +
			 w3 * *pT3 +
			 w4 * *pT4 +
			 w5 * *pT5 +
			 w6 * *pT6 +
			 w7 * *pT7
			 + pLayer->_cst28 * (
			 vT0 * *pL6d0s1 +
			 vT1 * *pL6d1s1 +
			 vT2 * *pL6d2s1 +
			 vT3 * *pL6d3s1 +
			 vT4 * *pL6d4s1 +
			 vT5 * *pL6d5s1 +
			 vT6 * *pL6d6s1 +
			 vT7 * *pL6d7s1 ));
	}


	return array;
}

double* CLevel6::update_noO_noDepth_gaussian(double* array, unsigned int direction, CLevel6* pLayer)
{
	double* p = array;

	double w0,w1,w2,w3,w4,w5,w6,w7;
	double v0,v1,v2,v3,v4,v5,v6,v7;

	unsigned int offset = pLayer->_nb_dir - direction;

	w0 = pLayer->_dg[offset];
	w1 = pLayer->_dg[offset+1];
	w2 = pLayer->_dg[offset+2];
	w3 = pLayer->_dg[offset+3];
	w4 = pLayer->_dg[offset+4];
	w5 = pLayer->_dg[offset+5];
	w6 = pLayer->_dg[offset+6];
	w7 = pLayer->_dg[offset+7];

	v0 = pLayer->_V[offset];
	v1 = pLayer->_V[offset+1];
	v2 = pLayer->_V[offset+2];
	v3 = pLayer->_V[offset+3];
	v4 = pLayer->_V[offset+4];
	v5 = pLayer->_V[offset+5];
	v6 = pLayer->_V[offset+6];
	v7 = pLayer->_V[offset+7];

	double* pT0 = pLayer->_baseV[0];
	double* pT1 = pLayer->_baseV[1];
	double* pT2 = pLayer->_baseV[2];
	double* pT3 = pLayer->_baseV[3];
	double* pT4 = pLayer->_baseV[4];
	double* pT5 = pLayer->_baseV[5];
	double* pT6 = pLayer->_baseV[6];
	double* pT7 = pLayer->_baseV[7];

	double* pL5d0 = pLayer->_baseRL5[0];
	double* pL5d1 = pLayer->_baseRL5[1];
	double* pL5d2 = pLayer->_baseRL5[2];
	double* pL5d3 = pLayer->_baseRL5[3];
	double* pL5d4 = pLayer->_baseRL5[4];
	double* pL5d5 = pLayer->_baseRL5[5];
	double* pL5d6 = pLayer->_baseRL5[6];
	double* pL5d7 = pLayer->_baseRL5[7];

	for(unsigned int i = 0; i < pLayer->_vSize * pLayer->_hSize; i++,p++,pT0++,pT1++,pT2++,pT3++,pT4++,pT5++,pT6++,pT7++,
														   pL5d0++,pL5d1++,pL5d2++,pL5d3++,pL5d4++,pL5d5++,pL5d6++,pL5d7++)
	{
		*p = pLayer->_cst27* *p + (1 - *p) * pLayer->_cst29 *
			(
			 v0 * *pL5d0 +
			 v1 * *pL5d1 +
			 v2 * *pL5d2 +
			 v3 * *pL5d3 +
			 v4 * *pL5d4 +
			 v5 * *pL5d5 +
			 v6 * *pL5d6 +
			 v7 * *pL5d7
			 ) +
			 (_cst30 + *p)	//No need to care about - sign,  dt or A9, since already weighted in the m_P kernel earlier
			 *
			 (w0 * *pT0 +
			 w1 * *pT1 +
			 w2 * *pT2 +
			 w3 * *pT3 +
			 w4 * *pT4 +
			 w5 * *pT5 +
			 w6 * *pT6 +
			 w7 * *pT7);
	}


	return array;
}

void CLevel6::print_kernels(const string& output_folder, bool recreate, string optName)
{
	KernelLayer::print_kernels(output_folder, recreate,optName);

	//1. Isotropic kernel
	string base = output_folder + "/" + _level_name;
	string endName = "IsotropicKernel";
	base += endName;
	base += optName;
	string endName2 = ".txt";
	base += endName2;

	_P.print_kernel(base, recreate);

	//2. V kernel
	base = output_folder + "/" + _level_name;
	base += "V.txt";
	ofstream vFile(base.c_str());
	vFile<<_V.size()<<endl;
	for(auto p = _V.begin() ; p != _V.end() ; p++)
		vFile<<*p<<endl;

	//3. V Tilde kernel
	base = output_folder + "/" + _level_name;
	base += "VTilde.txt";
	ofstream vtFile(base.c_str());
	vtFile<<_V_tilde.size()<<endl;
	for(auto p = _V_tilde.begin() ; p != _V_tilde.end() ; p++)
		vtFile<<*p<<endl;
}

void CLevel6::prune_kernel(Kernel<double>& k, double factor, double threshold)
{
	if(factor < threshold)
		k.remove();
	else
	{
		k.weight(factor);
		k.modulate_unit(factor);
	}

}

double CLevel6::top_down_weight(unsigned int x, unsigned int y, double a, double sigma, double attX, double attY)
{
	double num   = (x-attX)*(x-attX) + (y-attY)*(y-attY);
	double denom = sigma*sigma;
	return a * exp(-0.5 * num/denom);
}
}
