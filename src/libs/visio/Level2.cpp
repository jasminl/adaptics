#include "visio/Level2.h"
#include <fstream>

#pragma GCC diagnostic ignored "-Wsequence-point"	//Disable the sequence point warnings

using namespace std;

extern int		     g_countDownToSave;			    //Indicates whether to save current iteration or not

namespace visio
{
CLevel2::CLevel2(int iter_per_frame, double dt,
		    double A1,double B1,double C1,double A2,double K2,
			double A3,double B3,double C3,double K3,
			double A4,double B4,double C4,double K4,
			unsigned int verticalSize, unsigned int horizontalSize,
			unsigned int nbDirections,
			unsigned int totFrames,
			const string& output_folder,
			string name,
			bool record)
: CLayer(verticalSize, horizontalSize, totFrames, nbDirections,
		name, record, nullptr)
{
	//Pre-calculate parameters
	_cst1 = 1-dt*A1*B1;
	_cst2 = dt*A1*C1;
	_cst3 = -dt*A1;

	_cst4 = dt*A2;
	_cst5 = 1-dt*A2;
	_cst6 = -dt*A2*K2;

	_cst7 = 1 - dt*A3*B3;
	_cst8 = dt*A3*C3;
	_cst9 = -dt*A3*K3;

	_cst10 = 1 - dt*A4*B4;
	_cst11 = dt*A4*C4;
	_cst12 = -dt*A4*K4;

	//Create non-directional arrays
	create_array(_xOn, verticalSize, horizontalSize);
	create_array(_xOff, verticalSize, horizontalSize);
	create_array(_zOn, verticalSize, horizontalSize);
	create_array(_zOff, verticalSize, horizontalSize);
	create_array(_b, verticalSize, horizontalSize);

	//Set initial value of 1 to zOn and zOff
	set_value(_zOn, 1.0);
	set_value(_zOff, 1.0);

	//Create directional arrays
	create_array(_c0, verticalSize, horizontalSize);
	create_array(_c1, verticalSize, horizontalSize);
	create_array(_c2, verticalSize, horizontalSize);
	create_array(_c3, verticalSize, horizontalSize);
	create_array(_c4, verticalSize, horizontalSize);
	create_array(_c5, verticalSize, horizontalSize);
	create_array(_c6, verticalSize, horizontalSize);
	create_array(_c7, verticalSize, horizontalSize);

	create_array(_e0, verticalSize, horizontalSize);
	create_array(_e1, verticalSize, horizontalSize);
	create_array(_e2, verticalSize, horizontalSize);
	create_array(_e3, verticalSize, horizontalSize);
	create_array(_e4, verticalSize, horizontalSize);
	create_array(_e5, verticalSize, horizontalSize);
	create_array(_e6, verticalSize, horizontalSize);
	create_array(_e7, verticalSize, horizontalSize);

	if(keep_record())
		open_files(output_folder, iter_per_frame, dt);
}

CLevel2::~CLevel2()
{
	delete_array(_xOn);
	delete_array(_xOff);
	delete_array(_zOn);
	delete_array(_zOff);
	delete_array(_b);

	delete_array(_c0);
	delete_array(_c1);
	delete_array(_c2);
	delete_array(_c3);
	delete_array(_c4);
	delete_array(_c5);
	delete_array(_c6);
	delete_array(_c7);

	delete_array(_e0);
	delete_array(_e1);
	delete_array(_e2);
	delete_array(_e3);
	delete_array(_e4);
	delete_array(_e5);
	delete_array(_e6);
	delete_array(_e7);
}

double* CLevel2::array(const char* name)
{
	if(!strcmp(name,"XON"))		return _xOn;
	if(!strcmp(name,"XOFF"))	return _xOff;
	if(!strcmp(name,"ZON"))		return _zOn;
	if(!strcmp(name,"ZOFF"))	return _zOff;
	if(!strcmp(name,"B"))		return _b;

	if(!strcmp(name,"E0"))		return _e0;
	if(!strcmp(name,"E1"))		return _e1;
	if(!strcmp(name,"E2"))		return _e2;
	if(!strcmp(name,"E3"))		return _e3;
	if(!strcmp(name,"E4"))		return _e4;
	if(!strcmp(name,"E5"))		return _e5;
	if(!strcmp(name,"E6"))		return _e6;
	if(!strcmp(name,"E7"))		return _e7;

	return nullptr;
}

void CLevel2::open_files(const string& output_folder, int iter_per_frame, double dt)
{
	vector<string> ns(21);
	ns[0] = "XON.dat";ns[1] = "XOFF.dat";ns[2] = "ZON.dat";ns[3] = "ZOFF.dat";ns[4] = "B.dat";ns[5] = "C0.dat";ns[6] = "C1.dat";ns[7] = "C2.dat";
	ns[8] = "C3.dat";ns[9] = "C4.dat";ns[10] = "C5.dat";ns[11] = "C6.dat";ns[12] = "C7.dat";ns[13] = "E0.dat";ns[14] = "E1.dat";ns[15] = "E2.dat";
	ns[16] = "E3.dat";ns[17] = "E4.dat";ns[18] = "E5.dat";ns[19] = "E6.dat";ns[20] = "E7.dat";

	vector<ofstream*> fs(21);
	fs[0] = &_xOnStoreF;fs[1] = &_xOffStoreF;fs[2] = &_zOnStoreF;fs[3] = &_zOffStoreF;fs[4] = &_bStoreF;fs[5] = &_c0StoreF;fs[6] = &_c1StoreF;fs[7] = &_c2StoreF;
	fs[8] = &_c3StoreF;fs[9] = &_c4StoreF;fs[10] = &_c5StoreF;fs[11] = &_c6StoreF;fs[12] = &_c7StoreF;fs[13] = &_e0StoreF;fs[14] = &_e1StoreF;fs[15] = &_e2StoreF;
	fs[16] = &_e3StoreF;fs[17] = &_e4StoreF;fs[18] = &_e5StoreF;fs[19] = &_e6StoreF;fs[20] = &_e7StoreF;

	auto q = fs.begin();

	for(auto p = ns.begin(); p != ns.end() ; p++, q++)
	{
		string complete_name = output_folder + "/" + _level_name;
		complete_name += *p;
		(*q)->open(complete_name, ios::binary);

		if((*q)->bad())
			throw ios_base::failure("CLevel2::open_files: unable to open file" + complete_name);

		(*q)->write(reinterpret_cast<char*>(&_nb_frames), sizeof(unsigned int));
		(*q)->write(reinterpret_cast<char*>(/*&g_NbIterPerFrame*/&iter_per_frame), sizeof(unsigned int));
		(*q)->write(reinterpret_cast<char*>(&_vSize), sizeof(unsigned int));
		(*q)->write(reinterpret_cast<char*>(&_hSize), sizeof(unsigned int));
		(*q)->write(reinterpret_cast<char*>(/*&g_dt*/&dt), sizeof(double));
	}
}

void CLevel2::compute(const vector<int>& motionOn, const vector<int>& motionOff)
{
	unsigned int layer1_size = _hSize * _vSize;

	double* pOn   = _xOn;
	double* pOff  = _xOff;
	double* pZOn  = _zOn;
	double* pZOff = _zOff;
	double* pb    = _b;

	//Covariant pointers
	double* pc0   = _c0;
	double* pc1   = _c1 + _hSize;
	double* pc2   = _c2 + _hSize;
	double* pc3   = _c3 + _hSize+1;
	double* pc4   = _c4;
	double* pc5   = _c5;
	double* pc6   = _c6;
	double* pc7   = _c7;

	double* pb4   = _b + 1;
	double* pb1   = _b + _hSize;
	double* pb2   = _b + _hSize;
	double* pb3   = _b + _hSize+1;

	double* pe0   = _e0;
	double* pe1   = _e1 + _hSize;
	double* pe2   = _e2 + _hSize;
	double* pe3   = _e3 + _hSize+1;
	double* pe4   = _e4;
	double* pe5   = _e5;
	double* pe6   = _e6;
	double* pe7   = _e7;

	//Opponent pointers
	double* pc0D  = _c0;
	double* pc1D  = _c1+_hSize;
	double* pc2D  = _c2+_hSize;
	double* pc3D  = _c3+_hSize+1;
	double* pc4D  = _c4+1;
	double* pc5D  = _c5+1;
	double* pc6D  = _c6;
	double* pc7D  = _c7;

	static double xOnB, xOffB;
	static double cB[2] = {0, 0};

	auto pOnInput = motionOn.begin();
	auto pOffInput = motionOff.begin();

	unsigned int i;
	for(i = 0; i < layer1_size; i++, pOnInput++, pOffInput++, pOn++, pOff++, pb++, pZOn++, pZOff++)
	{
		//Calculate b first
		*pb = *pOn * *pZOn + *pOff * *pZOff;	//Note: here we could premultiply by m_cst8

		xOnB  = *pOn*_cst1 + _cst2 * *pOnInput + _cst3 * *pOn * *pOnInput;
		xOffB = *pOff*_cst1 + _cst2 * *pOffInput + _cst3 * *pOff * *pOffInput;

		//Update z
		*pZOn  = _cst5* *pZOn + _cst4 + _cst6 * *pOn * *pZOn;
		*pZOff  = _cst5* *pZOff + _cst4 + _cst6 * *pOff * *pZOff;

		//Update x
		*pOn = xOnB;
		*pOff= xOffB;
	}

	//Directions 0 and 4
	pb = _b;
	*pc4++ = _cst7 * *pc4 + _cst8 * *pb;
	*pe4++ = _cst10 * *pe4 + _cst11 * *pb;
	for(i = 0; i < layer1_size - 1; i++, pb++, pb4++, pc0++, pc0D++, pc4++, pc4D++, pe0++, pe4++)
	{
		if(*pc4D>0)
		{
			cB[0] = _cst7 * *pc0 + _cst8 * *pb + _cst9 * *pc4D;
			*pe0 = _cst10 * *pe0 + _cst11 * *pb + _cst12 * *pc4D;
		}
		else
		{
			cB[0] = _cst7 * *pc0 + _cst8 * *pb;
			*pe0 = _cst10 * *pe0 + _cst11 * *pb;
		}
		if(*pc0D>0)
		{
			cB[1] = _cst7 * *pc4 + _cst8 * *pb4 + _cst9 * *pc0D;
			*pe4 = _cst10 * *pe4 + _cst11 * *pb4 + _cst12 * *pc0D;
		}
		else
		{
			cB[1] = _cst7 * *pc4 + _cst8 * *pb4;
			*pe4 = _cst10 * *pe4 + _cst11 * *pb4;
		}

		*pc0 = cB[0];
		*pc4 = cB[1];

	}
 	*pc0 = _cst7 * *pc0 + _cst8 * *pb;
	*pe0 = _cst10 * *pe0 + _cst11 * *pb;

	//Directions 1 and 5
    pb = _b;
	*pc5++ = _cst7 * *pc5 + _cst8 * *pb;
	*pe5++ = _cst10 * *pe5 + _cst11 * *pb++;

	for(i = 0;i < layer1_size - _hSize - 1; i++, pb++, pb1++, pc1++, pc1D++, pc5++, pc5D++, pe1++, pe5++)
	{
		if(*pc5D>0)
		{
			cB[0] = _cst7 * *pc1 + _cst8 * *pb1 + _cst9 * *pc5D;
			*pe1 = _cst10 * *pe1 + _cst11 * *pb1 + _cst12 * *pc5D;
		}
		else
		{
			cB[0] = _cst7 * *pc1 + _cst8 * *pb1;
			*pe1 = _cst10 * *pe1 + _cst11 * *pb1;
		}

		if(*pc1D>0)
		{
			cB[1] = _cst7 * *pc5 + _cst8 * *pb + _cst9 * *pc1D;
			*pe5 = _cst10 * *pe5 + _cst11 * *pb + _cst12 * *pc1D;
		}
		else
		{
			cB[1] = _cst7 * *pc5 + _cst8 * *pb;
			*pe5 = _cst10 * *pe5 + _cst11 * *pb;
		}

		*pc1 = cB[0];
		*pc5 = cB[1];
	}
	*pc1++ = _cst7 * *pc1 + _cst8 * *pb1;
	*pe1++ = _cst10 * *pe1 + _cst11 * *pb1;
	memset(pc5, 0, _hSize * sizeof(double));	//Last row of direction 1 is undefined
	memset(_c1, 0, _hSize * sizeof(double));	//First row of direction 5 is undefined
	memset(pe5, 0, _hSize * sizeof(double));
	memset(_e1, 0, _hSize * sizeof(double));

	//Directions 2 and 6
	pb = _b;
	memset(_c2, 0, _hSize * sizeof(double));	//First row of direction 6 is undefined
	memset(_e2, 0, _hSize * sizeof(double));

	for(i = 0;i < layer1_size - _hSize; i++, pb++, pb2++, pc2++, pc2D++, pc6++, pc6D++, pe2++, pe6++)
	{
		if(*pc6D>0)
		{
			cB[0] = _cst7 * *pc2 + _cst8 * *pb2 + _cst9 * *pc6D;
			*pe2 = _cst10 * *pe2 + _cst11 * *pb2 + _cst12 * *pc6D;
		}
		else
		{
			cB[0] = _cst7 * *pc2 + _cst8 * *pb2 ;
			*pe2 = _cst10 * *pe2 + _cst11 * *pb2 ;
		}

		if(*pc2D>0)
		{
			cB[1] = _cst7 * *pc6 + _cst8 * *pb + _cst9 * *pc2D;
			*pe6 = _cst10 * *pe6 + _cst11 * *pb + _cst12 * *pc2D;
		}
		else
		{
			cB[1] = _cst7 * *pc6 + _cst8 * *pb  ;
			*pe6 = _cst10 * *pe6 + _cst11 * *pb ;
		}

		*pc2 = cB[0];
		*pc6 = cB[1];

	}

	memset(pc6, 0, _hSize * sizeof(double));	//Last row of direction 2 is undefined
	memset(pe6, 0, _hSize * sizeof(double));

	//Directions 3 and 7
	pb = _b;
	memset(_c3, 0, (_hSize+1) * sizeof(double));	//First row + one unit of direction 7 is undefined
	memset(_e3, 0, (_hSize+1) * sizeof(double));
	for(i = 0; i < layer1_size - _hSize - 1; i++, pb++, pb3++, pc3++, pc3D++, pc7++, pc7D++, pe7++, pe3++)
	{
		if(*pc7D>0)
		{
			cB[0] = _cst7 * *pc3 + _cst8 * *pb3 + _cst9 * *pc7D;
			*pe3 = _cst10 * *pe3 + _cst11 * *pb3 + _cst12 * *pc7D;

		}
		else
		{
			cB[0] = _cst7 * *pc3 + _cst8 * *pb3  ;
			*pe3 = _cst10 * *pe3 + _cst11 * *pb3  ;
		}
		if(*pc3D>0)
		{
			cB[1] = _cst7 * *pc7 + _cst8 * *pb + _cst9 * *pc3D;
			*pe7 = _cst10 * *pe7 + _cst11 * *pb + _cst12 * *pc3D;
		}
		else
		{
			cB[1] = _cst7 * *pc7 + _cst8 * *pb  ;
			*pe7 = _cst10 * *pe7 + _cst11 * *pb  ;
		}

		*pc3 = cB[0];
		*pc7 = cB[1];
	}
	memset(pc7, 0, _hSize * sizeof(double));	//Last row + one unit of direction 3 is undefined
	memset(pe7, 0, _hSize * sizeof(double));

	if(keep_record()) save();
}

void CLevel2::save()
{
	if(g_countDownToSave != 0) return;	//Save only if all steps are saved or only we're at the last one

	//Assumes file are already opened
	_xOnStoreF.write(reinterpret_cast<char*>(_xOn),sizeof(double)*_vSize*_hSize);
	_xOffStoreF.write(reinterpret_cast<char*>(_xOff),sizeof(double)*_vSize*_hSize);
	_zOnStoreF.write(reinterpret_cast<char*>(_zOn),sizeof(double)*_vSize*_hSize);
	_zOffStoreF.write(reinterpret_cast<char*>(_zOff),sizeof(double)*_vSize*_hSize);
	_bStoreF.write(reinterpret_cast<char*>(_b),sizeof(double)*_vSize*_hSize);
	_c0StoreF.write(reinterpret_cast<char*>(_c0),sizeof(double)*_vSize*_hSize);
	_c1StoreF.write(reinterpret_cast<char*>(_c1),sizeof(double)*_vSize*_hSize);
	_c2StoreF.write(reinterpret_cast<char*>(_c2),sizeof(double)*_vSize*_hSize);
	_c3StoreF.write(reinterpret_cast<char*>(_c3),sizeof(double)*_vSize*_hSize);
	_c4StoreF.write(reinterpret_cast<char*>(_c4),sizeof(double)*_vSize*_hSize);
	_c5StoreF.write(reinterpret_cast<char*>(_c5),sizeof(double)*_vSize*_hSize);
	_c6StoreF.write(reinterpret_cast<char*>(_c6),sizeof(double)*_vSize*_hSize);
	_c7StoreF.write(reinterpret_cast<char*>(_c7),sizeof(double)*_vSize*_hSize);
	_e0StoreF.write(reinterpret_cast<char*>(_e0),sizeof(double)*_vSize*_hSize);
	_e1StoreF.write(reinterpret_cast<char*>(_e1),sizeof(double)*_vSize*_hSize);
	_e2StoreF.write(reinterpret_cast<char*>(_e2),sizeof(double)*_vSize*_hSize);
	_e3StoreF.write(reinterpret_cast<char*>(_e3),sizeof(double)*_vSize*_hSize);
	_e4StoreF.write(reinterpret_cast<char*>(_e4),sizeof(double)*_vSize*_hSize);
	_e5StoreF.write(reinterpret_cast<char*>(_e5),sizeof(double)*_vSize*_hSize);
	_e6StoreF.write(reinterpret_cast<char*>(_e6),sizeof(double)*_vSize*_hSize);
	_e7StoreF.write(reinterpret_cast<char*>(_e7),sizeof(double)*_vSize*_hSize);
}
}
