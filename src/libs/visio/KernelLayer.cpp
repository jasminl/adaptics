#include "visio/KernelLayer.h"

using namespace std;

extern int		     g_countDownToSave;			    //Indicates whether to save current iteration or not

namespace visio
{
KernelLayer::KernelLayer(unsigned int VSize,unsigned int HSize, unsigned int totFrames, unsigned int nbDirections,
		   double sigma1Scale1, double sigma2Scale1, double sigma1Scale2, double sigma2Scale2,
		   string name, const string& output_folder, int iter_per_frame, double dt,
		   bool record, Buffer* buffer, double clipThreshold)
: CLayer(VSize, HSize, totFrames, nbDirections, name, record, buffer), _clip_thres(clipThreshold)
{
	if(sigma1Scale1 != NOT_ASSIGNED && sigma2Scale1 != NOT_ASSIGNED)
	{
		_d0s1K = Kernel<double>(0,sigma1Scale1,sigma2Scale1,clipThreshold);
		_d1s1K = Kernel<double>(PI/4,sigma1Scale1,sigma2Scale1,clipThreshold);
		_d2s1K = Kernel<double>(PI/2,sigma1Scale1,sigma2Scale1,clipThreshold);
		_d3s1K = Kernel<double>(3*PI/4,sigma1Scale1,sigma2Scale1,clipThreshold);
	}

	if(sigma1Scale2 != NOT_ASSIGNED && sigma2Scale2 != NOT_ASSIGNED)
	{
		_d0s2K = Kernel<double>(0,sigma1Scale2,sigma2Scale2,clipThreshold);
		_d1s2K = Kernel<double>(PI/4,sigma1Scale2,sigma2Scale2,clipThreshold);
		_d2s2K = Kernel<double>(PI/2,sigma1Scale2,sigma2Scale2,clipThreshold);
		_d3s2K = Kernel<double>(3*PI/4,sigma1Scale2,sigma2Scale2,clipThreshold);
	}

	_k = vector< Kernel<double>* > (8);
	_k[0] = &_d0s1K;
	_k[1] = &_d1s1K;
	_k[2] = &_d2s1K;
	_k[3] = &_d3s1K;
	_k[4] = &_d0s2K;
	_k[5] = &_d1s2K;
	_k[6] = &_d2s2K;
	_k[7] = &_d3s2K;

	//Create activation arrays
	CLayer::create_array(_d0s1,VSize,HSize);
	CLayer::create_array(_d1s1,VSize,HSize);
	CLayer::create_array(_d2s1,VSize,HSize);
	CLayer::create_array(_d3s1,VSize,HSize);
	CLayer::create_array(_d4s1,VSize,HSize);
	CLayer::create_array(_d5s1,VSize,HSize);
	CLayer::create_array(_d6s1,VSize,HSize);
	CLayer::create_array(_d7s1,VSize,HSize);

	CLayer::create_array(_d0s2,VSize,HSize);
	CLayer::create_array(_d1s2,VSize,HSize);
	CLayer::create_array(_d2s2,VSize,HSize);
	CLayer::create_array(_d3s2,VSize,HSize);
	CLayer::create_array(_d4s2,VSize,HSize);
	CLayer::create_array(_d5s2,VSize,HSize);
	CLayer::create_array(_d6s2,VSize,HSize);
	CLayer::create_array(_d7s2,VSize,HSize);

	_a = vector<double*> (16);
	_a[0] = _d0s1;
	_a[1] = _d1s1;
	_a[2] = _d2s1;
	_a[3] = _d3s1;
	_a[4] = _d4s1;
	_a[5] = _d5s1;
	_a[6] = _d6s1;
	_a[7] = _d7s1;
	_a[8] = _d0s2;
	_a[9] = _d1s2;
	_a[10] = _d2s2;
	_a[11] = _d3s2;
	_a[12] = _d4s2;
	_a[13] = _d5s2;
	_a[14] = _d6s2;
	_a[15] = _d7s2;

#if NUM_THREADS == 4 || NUM_THREADS == 16

	//Here, due to FFT buffers which can't be redundant across opposite direction diagonal kernels, create extra ones
	if(sigma1Scale1 != NOT_ASSIGNED && sigma2Scale1 != NOT_ASSIGNED)
	{
		_d5s1K = _d1s1K;
		_d7s1K = _d3s1K;
	}

	if(sigma1Scale2 != NOT_ASSIGNED && sigma2Scale2 != NOT_ASSIGNED)
	{
		_d5s2K = _d1s2K;
		_d7s2K = _d3s2K;
	}

	_k.insert(_k.end(),4,NULL);
	_k[8] = &_d5s1K;
	_k[9] = &_d7s1K;
	_k[10]= &_d5s2K;
	_k[11]= &_d7s2K;
/*
	m_k = vector< kernel<double>* > (12);
	m_k[0] = &m_d0s1K;
	m_k[1] = &m_d1s1K;
	m_k[2] = &m_d2s1K;
	m_k[3] = &m_d3s1K;
	m_k[4] = &m_d0s2K;
	m_k[5] = &m_d1s2K;
	m_k[6] = &m_d2s2K;
	m_k[7] = &m_d3s2K;
	m_k[8] = &m_d5s1K;
	m_k[9] = &m_d7s1K;
	m_k[10]= &m_d5s2K;
	m_k[11]= &m_d7s2K;
*/
#endif
	if(keep_record())
		open_files(output_folder, iter_per_frame, dt);
}

KernelLayer::~KernelLayer()
{
	//Delete activation values
	delete[] _d0s1;
	delete[] _d1s1;
	delete[] _d2s1;
	delete[] _d3s1;
	delete[] _d4s1;
	delete[] _d5s1;
	delete[] _d6s1;
	delete[] _d7s1;

	delete[] _d0s2;
	delete[] _d1s2;
	delete[] _d2s2;
	delete[] _d3s2;
	delete[] _d4s2;
	delete[] _d5s2;
	delete[] _d6s2;
	delete[] _d7s2;

	_d0s1 = _d1s1 = _d2s1 = _d3s1 = _d4s1 = _d5s1 = _d6s1 = _d7s1 = NULL;
	_d0s2 = _d1s2 = _d2s2 = _d3s2 = _d4s2 = _d5s2 = _d6s2 = _d7s2 = NULL;

	if(keep_record())
	{
		close_files();
	}

}


void KernelLayer::open_files(const std::string& output_folder, int iter_per_frame, double dt)
{
	vector<string> ns(16);
	ns[0] = "D0S1.dat";ns[1] = "D1S1.dat";ns[2] = "D2S1.dat";ns[3] = "D3S1.dat";ns[4] = "D4S1.dat";ns[5] = "D5S1.dat";ns[6] = "D6S1.dat";ns[7] = "D7S1.dat";
	ns[8] = "D0S2.dat";ns[9] = "D1S2.dat";ns[10] = "D2S2.dat";ns[11] = "D3S2.dat";ns[12] = "D4S2.dat";ns[13] = "D5S2.dat";ns[14] = "D6S2.dat";ns[15] = "D7S2.dat";

	vector<ofstream*> fs(16);
	fs[0] = &_d0s1F;fs[1] = &_d1s1F;fs[2] = &_d2s1F;fs[3] = &_d3s1F;fs[4] = &_d4s1F;fs[5] = &_d5s1F;fs[6] = &_d6s1F;fs[7] = &_d7s1F;
	fs[8] = &_d0s2F;fs[9] = &_d1s2F;fs[10] = &_d2s2F;fs[11] = &_d3s2F;fs[12] = &_d4s2F;fs[13] = &_d5s2F;fs[14] = &_d6s2F;fs[15] = &_d7s2F;

	auto q = fs.begin();

	for(auto p = ns.begin(); p != ns.end() ; p++, q++)
	{
		string completeName = output_folder + "/" + _level_name;
		completeName += *p;
		(*q)->open(completeName, ios::binary);

		(*q)->write(reinterpret_cast<char*>(&_nb_frames),sizeof(unsigned int));
		(*q)->write(reinterpret_cast<char*>(/*&g_NbIterPerFrame*/&iter_per_frame),sizeof(unsigned int));
		(*q)->write(reinterpret_cast<char*>(&_vSize),sizeof(unsigned int));
		(*q)->write(reinterpret_cast<char*>(&_hSize),sizeof(unsigned int));
		(*q)->write(reinterpret_cast<char*>(/*&g_dt*/&dt),sizeof(double));

		if((*q)->bad())
		{
			(*q)->close();
			return;
		}
	}
}

void KernelLayer::close_files()
{
	if(_d0s1F.is_open()) _d0s1F.close();
	if(_d1s1F.is_open()) _d1s1F.close();
	if(_d2s1F.is_open()) _d2s1F.close();
	if(_d3s1F.is_open()) _d3s1F.close();
	if(_d4s1F.is_open()) _d4s1F.close();
	if(_d5s1F.is_open()) _d5s1F.close();
	if(_d6s1F.is_open()) _d6s1F.close();
	if(_d7s1F.is_open()) _d7s1F.close();
	if(_d0s2F.is_open()) _d0s2F.close();
	if(_d1s2F.is_open()) _d1s2F.close();
	if(_d2s2F.is_open()) _d2s2F.close();
	if(_d3s2F.is_open()) _d3s2F.close();
	if(_d4s2F.is_open()) _d4s2F.close();
	if(_d5s2F.is_open()) _d5s2F.close();
	if(_d6s2F.is_open()) _d6s2F.close();
	if(_d7s2F.is_open()) _d7s2F.close();
}


void KernelLayer::time_kernel(unsigned int vsize, unsigned int hsize,Kernel<double>* k,const char* kname)
{
	string fname = _level_name;

	if(k != NULL)
	{
		fname.append(kname);
	}

	fname.append("ConvTime.txt");

	ofstream ofile(fname.c_str());

	//Create image and fill with random values
	double* randImage = new double[vsize*hsize];
	unsigned int i,j;
	for(i=0;i<vsize*hsize;i++)
	{
		randImage[i] =(double) rand()/RAND_MAX;
	}

	i=0;

	double* res = new double[vsize*hsize];
	double* buf = new double[vsize*hsize];

	clock_t start,end;

	const unsigned int nbTrials = 50;	//Do 10 trials of each and average
	unsigned int count=0,stdev=0;

	if(k != NULL)
	{
		//Test specific kernel
		k->svd();

		count = 0;
		stdev = 0;

		for(j=0;j<nbTrials;j++)
		{
			//Create new image
			for(unsigned int ki=0;ki<vsize*hsize;ki++)
			{
				randImage[ki] =(double) rand()/RAND_MAX;
			}

			start = clock();

			k->conv_sep(randImage,vsize,hsize,res,buf);	//Try SVD first

			end = clock();

			count += end-start;
			stdev += (end-start)*(end-start);
		}

		ofile<<"Kernel with SVD takes "<<(double)count/nbTrials/CLOCKS_PER_SEC<<" ("<<sqrt((double)stdev/nbTrials/CLOCKS_PER_SEC)<<") milliseconds, ";

		k->set_dft(vsize + k->vsize(), hsize + k->hsize());

		count = 0;
		stdev = 0;

		for(j=0;j<nbTrials;j++)
		{
			//Create new image
			for(unsigned int ki=0;ki<vsize*hsize;ki++)
			{
				randImage[ki] =(double) rand()/RAND_MAX;
			}

			start = clock();

			k->conv_dft(randImage,res,vsize,hsize);		//Try FFT second

			end = clock();

			count += end-start;
			stdev += (end-start)*(end-start);
		}

		ofile<<", with FFT takes "<<(double)count/nbTrials/CLOCKS_PER_SEC<<" ("<<sqrt((double)stdev/nbTrials/CLOCKS_PER_SEC)<<") milliseconds\t";

	}
	else
	{
		ofile<<"K\t"<<"SVD\t"<<"FFT\t";

		//Test all usual kernels
		for(vector< Kernel<double>* >::iterator p = _k.begin(); p != _k.end() ; p++,i++)
		{
			if((*p)->is_used())
			{
				if(!(i%2))
				{
					//Even direction (0,2,4,6): use SVD or FFT

					(*p)->svd();

					count = 0;
					stdev = 0;

					for(j=0;j<nbTrials;j++)
					{
						//Create new image
						for(unsigned int ki=0;ki<vsize*hsize;ki++)
						{
							randImage[ki] =(double) rand()/RAND_MAX;
						}

						start = clock();

						(*p)->conv_sep(randImage,vsize,hsize,res,buf);	//Try SVD first

						end = clock();

						count += end-start;
						stdev += (end-start)*(end-start);
					}

					ofile<<i<<"\t"<<(double)count/nbTrials/CLOCKS_PER_SEC<<" ("<<sqrt((double)stdev/nbTrials/CLOCKS_PER_SEC)<<")\t";

					(*p)->set_dft(vsize + (*p)->vsize(), hsize + (*p)->hsize());

					count = 0;

					for(j=0;j<nbTrials;j++)
					{
						//Create new image
						for(unsigned int ki=0;ki<vsize*hsize;ki++)
						{
							randImage[ki] =(double) rand()/RAND_MAX;
						}

						start = clock();

						(*p)->conv_dft(randImage,res,vsize,hsize);		//Try FFT second

						end = clock();

						count += end-start;
					}

					ofile<<(double)count/nbTrials/CLOCKS_PER_SEC<<" ("<<sqrt((double)stdev/nbTrials/CLOCKS_PER_SEC)<<")\t";

				}
				else
				{
					//Odd direction (1,3,5,7): use FFT or cudaFFT

					count = 0;

					(*p)->set_dft(vsize + (*p)->vsize(), hsize + (*p)->hsize());

					for(j=0;j<nbTrials;j++)
					{
						//Create new image
						for(unsigned int ki=0;ki<vsize*hsize;ki++)
						{
							randImage[ki] =(double) rand()/RAND_MAX;
						}

						start = clock();

						(*p)->conv_dft(randImage,res,vsize,hsize);

						end = clock();

						count += end-start;
					}

					ofile<<i<<"\t n.a.\t"<<(double)count/nbTrials/CLOCKS_PER_SEC<<" ("<<sqrt((double)stdev/nbTrials/CLOCKS_PER_SEC)<<")\t";
				}
			}
		}
	}

	ofile<<endl<<endl;

	ofile<<"m_k[0] = &m_d0s1K;\nm_k[1] = &m_d1s1K;\nm_k[2] = &m_d2s1K;\nm_k[3] = &m_d3s1K;\nm_k[4] = &m_d0s2K;\nm_k[5] = &m_d1s2K;\nm_k[6] = &m_d2s2K;\nm_k[7] = &m_d3s2K;\nm_k[8] = &m_d5s1K;\nm_k[9] = &m_d7s1K;\nm_k[10]= &m_d5s2K;\nm_k[11]= &m_d7s2K;";
	ofile.close();

	//Clear
	delete[] randImage;
	delete[] res;
	delete[] buf;
}


void KernelLayer::print_kernels(const string& output_folder, int recreate, const string optName)
{
	vector< Kernel<double>* > ks(16);
	ks[0] = &_d0s1K;ks[1] = &_d1s1K;ks[2]  = &_d2s1K;ks[3] = &_d3s1K;ks[4] = &_d0s1K;ks[5] = &_d1s1K;ks[6] = &_d2s1K;ks[7] = &_d3s1K;
	ks[8] = &_d0s2K;ks[9] = &_d1s2K;ks[10] = &_d2s2K;ks[11] = &_d3s2K;ks[12] = &_d0s2K;ks[13] = &_d1s2K;ks[14] = &_d2s2K;ks[15] = &_d3s2K;

	vector<string> ns(16);
	ns[0] = "D0S1.txt";ns[1] = "D1S1.txt";ns[2] = "D2S1.txt";ns[3] = "D3S1.txt";ns[4] = "D4S1.txt";ns[5] = "D5S1.txt";ns[6] = "D6S1.txt";ns[7] = "D7S1.txt";
	ns[8] = "D0S2.txt";ns[9] = "D1S2.txt";ns[10] = "D2S2.txt";ns[11] = "D3S2.txt";ns[12] = "D4S2.txt";ns[13] = "D5S2.txt";ns[14] = "D6S2.txt";ns[15] = "D7S2.txt";

	auto n = ns.begin();

	for(auto p = ks.begin() ; p != ks.end() ; p++, n++)
	{
		string base = output_folder + "/" + _level_name;
		string endName = "kernel";
		base += optName;	//Here will add an extra name tag only if optName is not default value ""
		base += endName;
		base += *n;
		(*p)->print_kernel(base, recreate);
	}
}

double* KernelLayer::array(const char* name)
{
	if(!strcmp(name,"D0S1"))		return _d0s1;
	if(!strcmp(name,"D1S1"))		return _d1s1;
	if(!strcmp(name,"D2S1"))		return _d2s1;
	if(!strcmp(name,"D3S1"))		return _d3s1;
	if(!strcmp(name,"D4S1"))		return _d4s1;
	if(!strcmp(name,"D5S1"))		return _d5s1;
	if(!strcmp(name,"D6S1"))		return _d6s1;
	if(!strcmp(name,"D7S1"))		return _d7s1;

	if(!strcmp(name,"D0S2"))		return _d0s2;
	if(!strcmp(name,"D1S2"))		return _d1s2;
	if(!strcmp(name,"D2S2"))		return _d2s2;
	if(!strcmp(name,"D3S2"))		return _d3s2;
	if(!strcmp(name,"D4S2"))		return _d4s2;
	if(!strcmp(name,"D5S2"))		return _d5s2;
	if(!strcmp(name,"D6S2"))		return _d6s2;
	if(!strcmp(name,"D7S2"))		return _d7s2;

	return NULL;
}

void KernelLayer::weight_kernels(double w)
{
	_d0s1K.weight(w);
	_d1s1K.weight(w);
	_d2s1K.weight(w);
	_d3s1K.weight(w);

	_d0s2K.weight(w);
	_d1s2K.weight(w);
	_d2s2K.weight(w);
	_d3s2K.weight(w);

}

void KernelLayer::weight_extern_kernel(Kernel<double>& k,double w)
{
	k.weight(w);
}

void KernelLayer::set_intrinsic_weight(double w)
{
	_d0s1K.set_intrinsic_weight(w);
	_d1s1K.set_intrinsic_weight(w);
	_d2s1K.set_intrinsic_weight(w);
	_d3s1K.set_intrinsic_weight(w);

	_d0s2K.set_intrinsic_weight(w);
	_d1s2K.set_intrinsic_weight(w);
	_d2s2K.set_intrinsic_weight(w);
	_d3s2K.set_intrinsic_weight(w);
}

void KernelLayer::set_intrinsic_weight_to(Kernel<double>& k, double w)
{
	k.set_intrinsic_weight(w);
}

void KernelLayer::delete_array(double** array,unsigned int nbrows)
{
	if(array == nullptr)	//Array was not assigned
		return;

	for(unsigned int i = 0; i < nbrows; i++)
		delete[] array[i];
	delete[] array;
}

void KernelLayer::create_array(double**& array, unsigned int nbrows, unsigned int nbcols)
{
	array = new double*[nbrows];
    for(unsigned int i=0;i<nbrows;i++)
	{
		array[i] = new double[nbcols];
		memset(array[i],0,sizeof(double)*nbcols);
	}
}

void KernelLayer::save()
{
	if(g_countDownToSave != 0 ) return;	//Save only if all steps are saved or only we're at the last one

	//Assumes file are already opened
	_d0s1F.write(reinterpret_cast<char*>(_d0s1),sizeof(double)*_vSize*_hSize);
	_d1s1F.write(reinterpret_cast<char*>(_d1s1),sizeof(double)*_vSize*_hSize);
	_d2s1F.write(reinterpret_cast<char*>(_d2s1),sizeof(double)*_vSize*_hSize);
	_d3s1F.write(reinterpret_cast<char*>(_d3s1),sizeof(double)*_vSize*_hSize);
	_d4s1F.write(reinterpret_cast<char*>(_d4s1),sizeof(double)*_vSize*_hSize);
	_d5s1F.write(reinterpret_cast<char*>(_d5s1),sizeof(double)*_vSize*_hSize);
	_d6s1F.write(reinterpret_cast<char*>(_d6s1),sizeof(double)*_vSize*_hSize);
	_d7s1F.write(reinterpret_cast<char*>(_d7s1),sizeof(double)*_vSize*_hSize);

	_d0s2F.write(reinterpret_cast<char*>(_d0s2),sizeof(double)*_vSize*_hSize);
	_d1s2F.write(reinterpret_cast<char*>(_d1s2),sizeof(double)*_vSize*_hSize);
	_d2s2F.write(reinterpret_cast<char*>(_d2s2),sizeof(double)*_vSize*_hSize);
	_d3s2F.write(reinterpret_cast<char*>(_d3s2),sizeof(double)*_vSize*_hSize);
	_d4s2F.write(reinterpret_cast<char*>(_d4s2),sizeof(double)*_vSize*_hSize);
	_d5s2F.write(reinterpret_cast<char*>(_d5s2),sizeof(double)*_vSize*_hSize);
	_d6s2F.write(reinterpret_cast<char*>(_d6s2),sizeof(double)*_vSize*_hSize);
	_d7s2F.write(reinterpret_cast<char*>(_d7s2),sizeof(double)*_vSize*_hSize);
}

double* KernelLayer::rectify(double* array, double* output,unsigned int mvsize, unsigned int mhsize)
{
	unsigned int sz = mvsize*mhsize;

	memcpy(output,array,sizeof(double)*sz);

	double* p = output;

	for(unsigned int i=0;i<sz;i++,p++)
	{
		if(*p<0) *p=0;
	}


	return output;
}

double* KernelLayer::rectify(double* array, double* output,unsigned int totSize)
{
	memcpy(output,array,sizeof(double) * totSize);

	double* p = output;

	for(unsigned int i=0;i<totSize;i++,p++)
	{
		if(*p<0) *p=0;
	}

	return output;
}

double* KernelLayer::rectify_square(double* array, double* output, unsigned int mvsize, unsigned int mhsize)
{
	unsigned int sz = mvsize*mhsize;
	memcpy(output,array,sizeof(double)*sz);

	double* p = output;

	for(unsigned int i=0;i<sz;i++,p++)
	{
		if(*p>0) *p *= *p;
		else
			*p = 0;
	}

	return output;
}
}
