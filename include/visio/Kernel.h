#pragma once

#include <iomanip>
#include <fstream>
#include <stdexcept>

namespace visio
{
/**
 *	Convolution kernels
 */
template <typename T>
class Kernel
{

public:
	//typedefs, etc.
	enum{NO_RECREATE_KER, RECREATE_KER};	//!<Whether to weight the kernel in print_kernel

public: //constructor

	/**
	 *	Default constructor
	 */
	Kernel<T>()
	{
		_kernel = nullptr;
		_vsize = _hsize = 0;
		_orientation = -1;
		_sigma1 = _sigma2 = -1;
		_intrinsic_weight = 1;
		_clipSize = 0;
		_mask1 = _mask2 = nullptr;
		_dft = _image = nullptr;
		_result = nullptr;
		_pad = nullptr;
	}

	/**
	 *	Constructor when size of kernel is unknown but clipThreshold is known.
	 *	\param orientation Orientation of the kernel.
	 *	\param sigma1 Spread along the principal axis.
	 *	\param sigma2 Spread along the minor axis.
	 *	\param clipThres Threshold below which values are set to zero.
	 */
	Kernel<T>(double orientation, double sigma1, double sigma2, double clipThres)
	:_orientation(orientation), _sigma1(sigma1), _sigma2(sigma2), _clipSize(clipThres)
	{
		determine_size(clipThres);
		make_kernel();

		_intrinsic_weight = 1;
		_mask1 = _mask2 = nullptr;
		_dft = _image = nullptr;
		_result = nullptr;

		_dftVSize = _dftHSize = 0;
		_pad = nullptr;

	}

	/**
	 *	Constructor when clip is unknown but size is known.
	 *	\param hsize Horizontal size of kernel.
	 *	\param vsize Vertical size of kernel.
	 *	\param orientation Orientation of kernel.
	 *	\param sigma1 Spread along the principal axis.
	 *	\param sigma2 Spread along the minor axis.
	 */
	Kernel<T>(unsigned int hsize, unsigned int vsize, double orientation, double sigma1, double sigma2)
	:_hsize(hsize), _vsize(vsize), _orientation(orientation), _sigma1(sigma1), _sigma2(sigma2)
	{
		if(is_even(hsize, vsize))
			throw std::runtime_error("kernel<T>::kernel: size must be odd.");

		make_kernel();

		_intrinsic_weight = 1;
		_clipSize = 0;

		_mask1 = _mask2 = nullptr;
		_dft = _image = nullptr;
		_result = _pad = nullptr;

		_dftVSize = _dftHSize = 0;

	}

	/**
	 *	Destructor
	 */
	virtual ~Kernel<T>()
	{
 		remove();
	}

protected: //Data members

	unsigned int _hsize;			//!< Horizontal size
	unsigned int _vsize;			//!< Vertical size
	double _orientation;			//!< Orientation (in rads)
	double _sigma1;				//!< Stdev of Gaussian
	double _sigma2;
	double _clipSize;				//!< Threshold value at which kernel is truncated

	T _intrinsic_weight;			//!< Weighting coefficient present in the kernel's equation

	T* _kernel;					//!< The kernel data

	T* _mask1;						//!< First SV mask (rows)
	T* _mask2;						//!< Second SV mask (cols)

	fftw_complex* _dft;			//!< FFT coefficients	for kernel
	fftw_complex* _image;			//!< FFT coefficients for convolved image
	fftw_plan _plan;				//!< FFT plan to generate kernel's FFT
	fftw_plan _upPlan;				//!< FFT plan to perform convolution
	fftw_plan _downPlan;			//!< FFT plan for inverse transform
	double* _result;				//!< Noncentered FFT output
	double* _pad;					//!< Padded array for input image
	unsigned int _dftVSize;		//!< Vertical size of padded array
	unsigned int _dftHSize;		//!< Horizontal size of padded array

public: //Member functions

	/**
		\return Vertical size of kernel.
	*/
	unsigned int vsize() const
	{
		return _vsize;
	}

	/**
		\return Horizontal size of kernel.
	*/
	unsigned int hsize() const
	{
		return _hsize;
	}

	/**
	 *	Frees the space of the kernel without destroying it.
	 */
	void remove()
	{
		delete[] _kernel;
		_kernel = NULL;

		//SVD related arrays
		delete[] _mask1;
		delete[] _mask2;
		_mask1 = _mask2 = NULL;

		//FFT related arrays
		if(_dft!=NULL)
		{
			fftw_free(_dft);
			_dft = NULL;
			fftw_free(_image);
			_image = NULL;
			fftw_destroy_plan(_plan);
			fftw_destroy_plan(_upPlan);
			fftw_destroy_plan(_downPlan);

			delete[] _result;
			delete[] _pad;
		}
	}

	/**
	 *	Sets the 2D center of the kernel to a particular value.
	 *	\param value Numeric value at which to set the kernel.
	 */
	void set_central(T value)
	{
		memset(_kernel, 0, sizeof(T) * _vsize * _hsize);
		_kernel[_hsize*(_vsize-1)/2 + (_hsize-1)/2] = value;

	}

	/**
	 *	\return A boolean flag indicating whether the kernel is used or not.
	 */
	bool is_used()
	{
		return (_kernel != nullptr);
	}

	/**
	 *	Indicates whether the sizes passed as argument are odd or even.
	 *	\param sz1 The first potential (vertical or horizontal) size to test.
	 *	\param sz2 The second potential (vertical or horizontal) size to test.
	 *	\return A Boolean flag indicating whether the proposed sizes satisfy the requirement that the kernel be odd.
	 */
	bool is_even(unsigned int sz1, unsigned int sz2) const
	{
		return (!(sz1%2) || !(sz2%2));
	}

	/**
	 *	Determines size of kernel (always a square).
	 *	\param thres the threshold below which weights are considered 0.
	 *	\return the size of one side of the kernel.
	 */
	unsigned int determine_size(T thres)
	{
		double sig = max(_sigma1, _sigma2);

		unsigned int loc = 0;

		T cT;

		do
		{
			cT = exp(-0.5 *  pow(loc / sig, 2));
			loc++;

		}while(cT > thres);
		loc--;

		return _hsize = _vsize = 2 * loc + 1;
	}

	/**
	 *	Additional scaling of the weight kernel.
	 *	\param w Weight to use to scale the kernel.
	 */
	void weight(T w)
	{
		T* p = _kernel;
		for(unsigned int i = 0; i < _hsize * _vsize; i++, p++)
		{
		//	*p++ = w * *p;//05/30/14 This was the original code and is replaced by the line below and putting p++ in the loop
			*p *= w;
		}
	}

	/**
	 *	To store the intrinsic weighting factor (the one that is specified in the kernel's equation).
	 *	\param w Intrinsic weight of the kernel.
	 */
	void set_intrinsic_weight(T w)
	{
		_intrinsic_weight = w;
	}

	/**
	 *	Regular 2D convolution.
	 *	\param matrix Input array to convolve.
	 *	\param mvsize Vertical size of the input.
	 *	\param mhsize Horizontal size of the input.
	 *	\param result Array containing the result of the convolution.
	 *	\return Pointer to the resulting array.
	 */
	T* convolve(T* matrix, unsigned int mvsize, unsigned mhsize, T* result)
	{
		unsigned int i, j, k, g;
		T* pRes = result;

		unsigned int hMargin = (_hsize - 1) / 2;
		unsigned int vMargin = (_vsize - 1) / 2;

		for(i = 0; i < mvsize; i++)
		{
			for(j = 0; j < mhsize; j++)
			{
				*pRes = 0;
				unsigned int dvBegin = min(vMargin,i);
				unsigned int dvEnd   = min(vMargin,mvsize-i-1);
                unsigned int dhBegin = min(hMargin,j);
				unsigned int dhEnd   = min(hMargin,mhsize-j-1);

                T* p = matrix + mhsize*(i-dvBegin) + j-dhBegin;
				T* q = &_kernel[(vMargin-dvBegin)*_hsize + hMargin-dhBegin];

				unsigned int vDo    = dvEnd+dvBegin+1;
				unsigned int hDo    = dhEnd+dhBegin+1;

				unsigned int hSkip  = _hsize - hDo;
				unsigned int mhSkip = mhsize - hDo;

				for(k=0;k<vDo;k++)
				{
					for(g=0;g<hDo;g++,p++,q++)
					{
						*pRes += *q * *p;

					}

					q+= hSkip;
					p+= mhSkip;
				}

				pRes++;
			}
		}

		return result;
	}

	/**
	 *	Convolution with additional rectification
	 *	\param matrix Input array to convolve.
	 *	\param mvsize Vertical size of the input.
	 *	\param mhsize Horizontal size of the input.
	 *	\param result Array containing the result of the convolution.
	 *	\return Pointer to the resulting array.
	 */
	T* convolve_rectified(T* matrix, unsigned int mvsize, unsigned mhsize, T* result)
	{
		unsigned int i,j,k,g;
		T* pRes = result;

		unsigned int hMargin = (_hsize-1)/2;
		unsigned int vMargin = (_vsize-1)/2;

		for(i=0;i<mvsize;i++)
		{
			for(j=0;j<mhsize;j++)
			{
				*pRes = 0;

				unsigned int dvBegin = min(vMargin,i);
				unsigned int dvEnd   = min(vMargin,mvsize-i-1);
                unsigned int dhBegin = min(hMargin,j);
				unsigned int dhEnd   = min(hMargin,mhsize-j-1);

                T* p = matrix + mhsize*(i-dvBegin) + j-dhBegin;
				T* q = &_kernel[(vMargin-dvBegin)*_hsize + hMargin-dhBegin];

				unsigned int vDo    = dvEnd+dvBegin+1;
				unsigned int hDo    = dhEnd+dhBegin+1;

				unsigned int hSkip  = _hsize - hDo;
				unsigned int mhSkip = mhsize - hDo;

				for(k=0;k<vDo;k++)
				{
					for(g=0;g<hDo;g++,p++,q++)
					{
						if( *p > 0 )	//This is the only difference with the non-rectified version
							*pRes += *q * *p;
					}

					q+= hSkip;
					p+= mhSkip;
				}

				pRes++;
			}
		}

		return result;
	}


	/**
	 *	Convolution with additional rectification and squaring
	 *	\param matrix Input array to convolve.
	 *	\param mvsize Vertical size of the input.
	 *	\param mhsize Horizontal size of the input.
	 *	\param result Array containing the result of the convolution.
	 *	\return Pointer to the resulting array.
	 */
	T* convolve_rectified_squared(T* matrix, unsigned int mvsize,unsigned mhsize,T* result)
	{
		unsigned int i,j,k,g;
		T* pRes = result;

		unsigned int hMargin = (_hsize-1)/2;
		unsigned int vMargin = (_vsize-1)/2;

		for(i=0;i<mvsize;i++)
		{
			for(j=0;j<mhsize;j++)
			{
				*pRes = 0;

				unsigned int dvBegin = min(vMargin,i);
				unsigned int dvEnd   = min(vMargin,mvsize-i-1);
                unsigned int dhBegin = min(hMargin,j);
				unsigned int dhEnd   = min(hMargin,mhsize-j-1);

                T* p = matrix + mhsize*(i-dvBegin) + j-dhBegin;
				T* q = &_kernel[(vMargin-dvBegin)*_hsize + hMargin-dhBegin];

				unsigned int vDo    = dvEnd+dvBegin+1;
				unsigned int hDo    = dhEnd+dhBegin+1;

				unsigned int hSkip  = _hsize - hDo;
				unsigned int mhSkip = mhsize - hDo;

				for(k=0;k<vDo;k++)
				{
					for(g=0;g<hDo;g++,p++,q++)
					{
						if( *p > 0 )	//This is the only difference with the previous version
							*pRes += *q * *p * *p;
					}

					q+= hSkip;
					p+= mhSkip;
				}

				pRes++;
			}
		}

		return result;
	}

	/**
	 *	Creates the kernel.
	 *	\param buffer If NULL, creates m_kernel, otherwise uses buffer to store weights.
	 *	\param weight This parameter is used to further manually weight the kernel (useful to weight by intrinsic weight when recreating).
	 */
	void make_kernel(T* buffer=NULL,T weight=1)
	{
		T* pker;

		if(buffer==NULL)
		{
			_kernel = new T[_hsize*_vsize];
			pker = _kernel;
		}
		else
		{
			pker = buffer;
		}

		double cx;
		double cy;

		if(!(_hsize%2))	cx = _hsize/2.0 - 0.5;			//Even size, shift center by half unit
		else				cx = (_hsize-1)/2.0 ;		//Odd size

		if(!(_vsize%2))	cy = _vsize/2.0 - 0.5;			//Even size, shift center by half unit
		else				cy = (_vsize-1)/2.0 ;		//Odd size

		double rot[2][2];
		rot[0][0] = cos(_orientation);
		rot[0][1] = -sin(_orientation);
		rot[1][0] = sin(_orientation);
		rot[1][1] = cos(_orientation);

		static double ni,nj;

		//Note: j represents x coordinate, i represents y coordinate
		for(unsigned int i=0;i<_vsize;i++)
		{
			for(unsigned int j=0;j<_hsize;j++)
			{
				nj = rot[0][0] * (j-cx) + rot[0][1] * (i-cy);
				ni = rot[1][0] * (j-cx) + rot[1][1] * (i-cy);

				pker[i*_hsize + j] =  weight * static_cast<T>(exp(-0.5 * ( pow(ni/_sigma2,2) + pow(nj/_sigma1,2)  )));
			}
		}

	}

	/**
	 *	To print out kernel to stdout.
	 *	\param precision Degree of numerical precision desired.
	 */
	void show_kernel(unsigned int precision) const
	{
		for(unsigned int i=0;i<_vsize;i++)
		{
			for(unsigned int j=0;j<_hsize;j++)
			{
				std::cout<<std::setprecision(precision)<<_kernel[i*_hsize + j]<<" ";
			}
			cout<<endl;
		}
	}

	/**
	 *	Outputs kernel to a file.
	 *	\param fname Name of file to which to print the kernel.
	 *	\param recreate Flag that indicates whether the kernel should be printed out without the intrinsic weight factoring.
	 */
	void print_kernel(string fname, int recreate)
	{
		if (_kernel == nullptr)
			return;		//Signifies the kernel is nonexistent

		std::ofstream ofile(fname);

		if(recreate == RECREATE_KER)
		{
			//Recreate kernel to print it without the weighting, except for intrinsic weight
			T* buffer = new T[_hsize*_vsize];

			make_kernel(buffer, _intrinsic_weight);
			for(unsigned int i = 0; i < _vsize; i++)
			{
				for(unsigned int j = 0; j < _hsize; j++)
					ofile<<buffer[i*_hsize + j]<<" ";
				ofile<<endl;
			}
			delete[] buffer;
		}
		else
		{
			for(unsigned int i = 0; i < _vsize; i++)
			{
				for(unsigned int j = 0; j < _hsize; j++)
					ofile<<_kernel[i*_hsize + j]<<" ";
				ofile<<endl;
			}
		}
	}

	/**
	 *	Weights a particular connection.
	 *	\param term Quantity to add to the modulated weight.
	 *	\param row Index of the row of the unit to modulate.
	 *	\param column Index of the columns of the unit to modulate.
	 */
	void modulate_unit(T term, int row=-1, int column=-1)
	{
		if(row==-1 || column==-1)
		{
			//These particular indexes mean that the center unit must be modulated
			row    = _vsize/2;
			column = _hsize/2;
		}

		_kernel[row*_hsize + column] += term;
	}

	/**
	 *	Performs svd of kernel (for faster separable convolution).
	 *	Uses the JAMA library (http://math.nist.gov/tnt/download.html)
	 */
	void svd()
	{

		if(_kernel==NULL)
		{
			//Either the kernel hasn't been assigned or it has been pruned, so no way to do svd
			return;
		}

		//Create TNT/JAMA objects
		TNT::Array2D<T> a(_vsize,_hsize);

		T* p = _kernel;
		T** q = a;

		unsigned int i,j;

		//Assign kernel values
		for(i=0;i<_vsize;i++)
		{
			for(j=0;j<_hsize;j++,p++)
				q[i][j] = *p;
		}

  		JAMA::SVD<T> b(a);
		TNT::Array2D<T> u,v;
		TNT::Array1D<T> s;
		b.getU(u);
		b.getV(v);
		b.getSingularValues(s);

		//Create m_mask1 and m_mask2
		_mask1 = new T[u.dim2()];
		_mask2 = new T[v.dim1()];

		T* pm;
		for(i=0,pm = _mask1,q = u;i<static_cast<unsigned>(u.dim2());i++,pm++) *pm = q[i][0] * s[0];
		for(i=0,pm = _mask2,q = v;i<static_cast<unsigned>(v.dim1());i++,pm++) *pm = q[i][0];

	}

	/**
	 *	DFT setup
	 *	\param mvise Vertical size of the typical input frame.
	 *	\param mhsize Horizontal size of the typical input frame.
	 */
	void set_dft(unsigned int mvsize, unsigned int mhsize);

	/**
	 *	To center the kernel such that middle is in location [0 0] for FFT convolutions.
	 *	\param ker Kernel to center.
	 *	\param kerVSize Vertical size of the kernel.
	 *	\param kerHSize Horizontal size of the kernel.
	 *	\param output Array containing the padded centered kernel.
	 *	\param vSize Vertical size of the array containing the padded centered kernel.
	 *	\param hSize Horizontal size of the array containing the padded centered kernel.
	 */
	void center_kernel(double* ker, unsigned int kerVSize, unsigned int kerHSize, double* output, unsigned int vSize, unsigned int hSize)
	{
		unsigned int hbound = (kerHSize-1)/2;
		unsigned int vbound = (kerVSize-1)/2;

		unsigned int i;
		double* p = output;
		double* q = ker;

		//1. Copy lower right corner
		for(i=0;i<vbound+1;i++)
		{
			memcpy(p,q+(kerVSize*i)+kerVSize*hbound+hbound,sizeof(double) * (hbound+1));
			p+=hSize;
		}

		//2. Copy upper right corner
		p = output+hSize*(vSize-vbound);
		for(i=0;i<vbound;i++)
		{
			memcpy(p,q+(kerVSize*i)+hbound,sizeof(double) * (hbound+1));
			p+= hSize;
		}

		//3. Copy lower left corner
		p = output + hSize-hbound;
		for(i=0;i<vbound+1;i++)
		{
			memcpy(p,q+(kerVSize*i)+kerVSize*hbound,sizeof(double) * hbound);
			p+= hSize;
		}

		//4. Copy upper left corner
		p = output + hSize*(vSize-vbound) + (hSize-hbound);
		for(i=0;i<vbound;i++)
		{
			memcpy(p,q+kerVSize*i,sizeof(double)*hbound);
			p+= hSize;
		}
	}

	/**
	 *	FFT convolution
	 *	\param image Input frame to convolve.
	 *	\param result Array in which to store the result.
	 *	\param vsize Vertical size of the input frame.
	 *	\param hsize Horizontal size of the input frame.
	 */
	void conv_dft(double* image, double* result,unsigned vsize, unsigned hsize)
	{
		unsigned int i,j,ij=0;

		double* pim = image;
		double* ppa = _pad;

		//TODO: remove some of these?
		memset(_pad,0,sizeof(double)*_dftVSize*_dftHSize);
		memset(_result,0,sizeof(double)*_dftVSize*_dftHSize);
		memset(_image,0,sizeof(double)*_dftVSize*(_dftHSize/2+1));

		//0. Copy into padded array
		for(i=0;i<vsize;i++)
		{
			memcpy(ppa,pim,sizeof(double)*hsize);
			pim  += hsize;
			ppa  += _dftHSize;
		}

		//1. Perform DFT
		fftw_execute(_upPlan);

		double buf;
		double scale = 1.0/(_dftVSize*_dftHSize);

		ij=0;
		for(i=0;i<_dftVSize;i++)
		{
			for(j=0;j<_dftHSize/2+1;j++)
			{
				buf = (_image[ij][0]*_dft[ij][0] - _image[ij][1]*_dft[ij][1])*scale;
				_image[ij][1] = (_image[ij][0]*_dft[ij][1] + _image[ij][1]*_dft[ij][0])*scale;
				_image[ij][0] = buf;

				ij++;
			}
		}

		//3. Inverse DFT
		fftw_execute(_downPlan);

		//4. Copy to result array
		double* pres = result;
		double* pinv = _result;

		for(i=0;i<vsize;i++)
		{
			memcpy(pres,pinv,sizeof(double)*hsize);
			pres += hsize;
			pinv += _dftHSize;
		}

	}

	/**
	 *	Separable convolution. Convolves using m_mask1 and m_mask2 for separable kernels.
	 *	\param matrix Input frame.
	 *	\param mvsize Vertical size of the input frame.
	 *	\param mhsize Horizontal size of the input frame.
	 *	\param result Array storing the result of the convolution.
	 *	\param buffer Array used as buffer as part of the computations.
	 *	\return Pointer to the array storing the result.
	 */
	T* conv_sep(T* matrix, unsigned int mvsize,unsigned mhsize,T* result,T* buffer)
	{
		T* pRes = buffer;
		unsigned int i,j,k;

		//Convolve 1D on rows, store in 'result'
		unsigned int hMargin = (_hsize-1)/2;
		unsigned int vMargin = (_vsize-1)/2;

		for(i=0;i<mvsize;i++)
		{
			for(j=0;j<mhsize;j++)
			{
				*pRes = 0;

                unsigned int dhBegin = min(hMargin,j);
				unsigned int dhEnd   = min(hMargin,mhsize-j-1);

                T* p = matrix + mhsize*i + j-dhBegin;
				T* q = &_mask2[hMargin-dhBegin];

				unsigned int hDo    = dhEnd+dhBegin+1;

				for(k=0;k<hDo;k++,p++,q++)
				{
					*pRes += *q * *p;
				}

				pRes++;
			}

		}

		pRes = result;

		//Convolve 1D on cols, store in 'result'
		for(i=0;i<mvsize;i++)
		{
			for(j=0;j<mhsize;j++)
			{
				*pRes = 0;

				unsigned int dvBegin = min(vMargin,i);
				unsigned int dvEnd   = min(vMargin,mvsize-i-1);

				T* p = buffer + mhsize*(i-dvBegin) + j;
				T* q = &_mask1[vMargin-dvBegin];

				unsigned int vDo    = dvEnd+dvBegin+1;
				unsigned int mhSkip = mhsize;//hDo;

				for(k=0;k<vDo;k++,q++,p+= mhSkip)
				{
					*pRes += *q * *p;
				}

				pRes++;
			}
		}

		return result;
	}

	/**
	 *	Convolves with rectification using m_mask1 and m_mask2 for separable kernels.
	 *	\param matrix Input frame.
	 *	\param mvsize Vertical size of the input frame.
	 *	\param mhsize Horizontal size of the input frame.
	 *	\param result Array storing the result of the convolution.
	 *	\param buffer Array used as buffer as part of the computations.
	 *	\return Pointer to the array storing the result.
	 */
	T* conv_sep_rectified(T* matrix, unsigned int mvsize,unsigned mhsize,T* result,T* buffer)
	{
		T* pRes = buffer;
		unsigned int i,j,k;

		//Convolve 1D on rows, store in 'result'
		unsigned int hMargin = (_hsize-1)/2;
		unsigned int vMargin = (_vsize-1)/2;

		for(i=0;i<mvsize;i++)
		{
			for(j=0;j<mhsize;j++)
			{
				*pRes = 0;

                unsigned int dhBegin = min(hMargin,j);
				unsigned int dhEnd   = min(hMargin,mhsize-j-1);

                T* p = matrix + mhsize*i + j-dhBegin;
				T* q = &_mask2[hMargin-dhBegin];

				unsigned int hDo    = dhEnd+dhBegin+1;

				for(k=0;k<hDo;k++,p++,q++)
				{
					if(*p>0)
					*pRes += *q * *p;
				}

				pRes++;
			}

		}

		pRes = result;

		//Convolve 1D on cols, store in 'result'
		for(i=0;i<mvsize;i++)
		{
			for(j=0;j<mhsize;j++)
			{
				*pRes = 0;

				unsigned int dvBegin = min(vMargin,i);
				unsigned int dvEnd   = min(vMargin,mvsize-i-1);

				T* p = buffer + mhsize*(i-dvBegin) + j;
				T* q = &_mask1[vMargin-dvBegin];

				unsigned int vDo    = dvEnd+dvBegin+1;
				unsigned int mhSkip = mhsize;//hDo;

				for(k=0;k<vDo;k++,q++,p+= mhSkip)
				{
					*pRes += *q * *p;
				}

				pRes++;
			}
		}

		return result;
	}

	/**
	 *	Convolves with rectification and squaring using m_mask1 and m_mask2 for separable kernels
	 *	\param matrix Input frame.
	 *	\param mvsize Vertical size of the input frame.
	 *	\param mhsize Horizontal size of the input frame.
	 *	\param result Array storing the result of the convolution.
	 *	\param buffer Array used as buffer as part of the computations.
	 *	\return Pointer to the array storing the result.
	 */
	T* conv_sep_rectified_squared(T* matrix, unsigned int mvsize,unsigned mhsize,T* result,T* buffer)
	{
		T* pRes = buffer;
		unsigned int i,j,k;

		//Convolve 1D on rows, store in 'result'
		unsigned int hMargin = (_hsize-1)/2;
		unsigned int vMargin = (_vsize-1)/2;

		for(i=0;i<mvsize;i++)
		{
			for(j=0;j<mhsize;j++)
			{
				*pRes = 0;

                unsigned int dhBegin = min(hMargin,j);
				unsigned int dhEnd   = min(hMargin,mhsize-j-1);

                T* p = matrix + mhsize*i + j-dhBegin;
				T* q = &_mask2[hMargin-dhBegin];

				unsigned int hDo    = dhEnd+dhBegin+1;

				for(k=0;k<hDo;k++,p++,q++)
				{
					if(*p>0)
					*pRes += *q * *p * *p;
				}

				pRes++;
			}

		}

		pRes = result;

		//Convolve 1D on cols, store in 'result'
		for(i=0;i<mvsize;i++)
		{
			for(j=0;j<mhsize;j++)
			{
				*pRes = 0;

				unsigned int dvBegin = min(vMargin,i);
				unsigned int dvEnd   = min(vMargin,mvsize-i-1);

				T* p = buffer + mhsize*(i-dvBegin) + j;
				T* q = &_mask1[vMargin-dvBegin];

				unsigned int vDo    = dvEnd+dvBegin+1;
				unsigned int mhSkip = mhsize;//hDo;

				for(k=0;k<vDo;k++,q++,p+= mhSkip)
				{
					*pRes += *q * *p ;
				}

				pRes++;
			}
		}

		return result;
	}



public:

	/**
	 *	Assigns the value of a kernel to another.
	 *	\param a Base kernel from which to create a new kernel.
	 */
	Kernel<T>& operator=(const Kernel<T>& a)
	{
		_orientation = a._orientation;
		_sigma1 = a._sigma1;
		_sigma2 = a._sigma2;
		_hsize = a._hsize;
		_vsize = a._vsize;

		_kernel = new T[_hsize*_vsize];

		memcpy(_kernel,a._kernel,sizeof(T)* _hsize * _vsize);

		_pad = nullptr;

		set_intrinsic_weight(a._intrinsic_weight);
		return *this;
	}
};

/**
 * Class for convolution kernels.
 */
template<typename T>
inline void Kernel<T>::set_dft(unsigned int mvsize, unsigned int mhsize)
{
	if (_kernel == nullptr)
	{
		//This signifies either that the kernel hasn't been allocated or that it has been pruned
		return;
	}
	//0. Set sizes for convDFT()
	_dftVSize = mvsize;
	_dftHSize = mhsize;
	//1. Make kernel plan and compute coefficients
	double* b = new double[mvsize * mhsize];
	_dft = (fftw_complex*) (fftw_malloc(
			sizeof(fftw_complex) * mvsize * (mhsize / 2 + 1)));
	_plan = fftw_plan_dft_r2c_2d(mvsize, mhsize, b, _dft, FFTW_MEASURE);
	memset(b, 0, sizeof(double) * mvsize * mhsize);
	center_kernel(_kernel, _vsize, _hsize, b, mvsize, mhsize);
	fftw_execute(_plan);
	delete[] b;
	_pad = new double[mvsize * mhsize];
	//2. Make image plan
	_image = (fftw_complex*) (fftw_malloc(
			sizeof(fftw_complex) * mvsize * (mhsize / 2 + 1)));
	_result = new double[mvsize * mhsize];
	_upPlan = fftw_plan_dft_r2c_2d(mvsize, mhsize, _pad, _image, FFTW_MEASURE);
	_downPlan = fftw_plan_dft_c2r_2d(mvsize, mhsize, _image, _result,
			FFTW_MEASURE);
	memset(_pad, 0, sizeof(double) * mvsize * mhsize);
}
}
