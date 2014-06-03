#pragma once

#include "Layer.h"
#include "Kernel.h"

namespace visio
{
/**
	Kernel layer base class
*/
class KernelLayer:public CLayer
{

public:
	/**
	 *	Constructor
	 */
	KernelLayer(unsigned int VSize,unsigned int HSize, unsigned int totFrames, unsigned int nbDirections,
		   double sigma1Scale1, double sigma2Scale1, double sigma1Scale2, double sigma2Scale2,
		   const std::string name, const std::string& output_folder, int iter_per_frame, double dt,
		   bool record = false, Buffer* buffer = nullptr,
		   double clipThreshold = 0.01);

	/**
	 *	Destructor
	 */
	virtual
	~KernelLayer();

	/** @name Access methods
		These methods are used to access cell activities directly.
	*/
	//@{
	double* d0s1(){return _d0s1;}
	double* d1s1(){return _d1s1;}
	double* d2s1(){return _d2s1;}
	double* d3s1(){return _d3s1;}
	double* d4s1(){return _d4s1;}
	double* d5s1(){return _d5s1;}
	double* d6s1(){return _d6s1;}
	double* d7s1(){return _d7s1;}
	double* d0s2(){return _d0s2;}
	double* d1s2(){return _d1s2;}
	double* d2s2(){return _d2s2;}
	double* d3s2(){return _d3s2;}
	double* d4s2(){return _d4s2;}
	double* d5s2(){return _d5s2;}
	double* d6s2(){return _d6s2;}
	double* d7s2(){return _d7s2;}
	//@}

	/**
	 * To print all kernels to files. If recreate is true, it recreates the weights in a external buffer (useful to
	 * display weights without weighting factor) before printing them. Optname is useful if you want to add a specific
	 * part to the generated filenames.
	 * \param recreate Indicates whether to include intrinsic weighting or not.
	 * \param optName Optional file name where to save the kernels.
	 */
	virtual
	void print_kernels(const std::string& output_folder, int recreate = Kernel<double>::NO_RECREATE_KER,
			const std::string optName= "");

	/**
	 *	Access activity array by name.
	 *	\param name Name of the array to access.
	 *	\return Pointer to the desired array.
	 */
	double* array(const char* name);

	/**
	 *	For pre-weighting kernels.
	 *	\param Weight to use as factor.
	 */
	void weight_kernels(double w);

	/**
	 *	For pre-weighting a kernel not part of the layer.
	 *	\param k External kernel to weight.
	 *	\param w Weight to use as factor.
	 */
	void weight_extern_kernel(Kernel<double>& k, double w);

	/**
	 *	Stores weighting factor to allow recovering kernel without/with it later on (for display).
	 *	\param w Weight to use as factor.
	 */
	void set_intrinsic_weight(double w);

	/**
	 *	Stores weighting factor to allow recovering of an external kernel without/with it later on (for display).
	 *	\param k External kernel to weight.
	 *	\param w Weight to use as factor.
	 */
	void set_intrinsic_weight_to(Kernel<double>& k, double w);

	/**
	 *	Frees the space associated with a 2D array.
	 *	\param array Array to delete.
	 *	\param nbrows The number of rows in the array.
	 */
	void delete_array(double** array, unsigned int nbrows);

	/**
	 *	Allocates space for a 2D array.
	 *	\param array Pointer to array to be allocated.
	 *	\param nbrows Number of rows.
	 *	\param nbcols Number of columns.
	 */
	void create_array(double**& array, unsigned int nbrows, unsigned int nbcols);

	/** @name Overloaded files for saving results
	*/
	//@{
	void save();

	/**
	 * Overloaded function to save activity.
	 */
	void open_files(const std::string& output_folder, int iter_per_frame, double dt);

	/**
	 * Overloaded method to save files
	 */
	void close_files();
	//@}

	/**
	 *	Measures time for a given kernel on an input image
	 *	\param vsize Vertical size of hypothetical input frame.
	 *	\param hsize Horizontal size of a hypothetical input frame.
	 *	\param k Optional external kernel to time.
	 *	\param kname Optional name of the kernel (to store in the results).
	 */
	void time_kernel(unsigned int vsize, unsigned int hsize, Kernel<double>* k = nullptr, const char* kname = nullptr);

	/** @name Rectification methods
	 *	Rectification methods with externally provided buffer output.
	 */
	//@{
	double* rectify(double* array, double* output,unsigned int mvsize, unsigned int mhsize);

	/**
	 * Basic half-wave rectification
	 */
	static
	double* rectify(double* array, double* output,unsigned int totSize);

	/**
	 * Half-wave rectification with squaring
	 */
	double* rectify_square(double* array, double* output, unsigned int mvsize, unsigned int mhsize);
	//@}

protected:	//Data members

	/** @name Kernels
	*/
	//@{
	Kernel<double> _d0s1K;		//!< Direction 0, scale 1
	Kernel<double> _d1s1K;		//!< Direction 1, scale 1
	Kernel<double> _d2s1K;		//!< Direction 2, scale 1
	Kernel<double> _d3s1K;		//!< Direction 3, scale 1

	Kernel<double> _d0s2K;		//!< Direction 0, scale 2
	Kernel<double> _d1s2K;		//!< Direction 1, scale 2
	Kernel<double> _d2s2K;		//!< Direction 2, scale 2
	Kernel<double> _d3s2K;		//!< Direction 3, scale 2
	//@}

	std::vector<double*> _a;			//!< Array of pointers to activation variables
	std::vector< Kernel<double>* > _k;	//!< Array of pointers to kernels

#if NUM_THREADS==4 || NUM_THREADS==16

	/** @name Supplementary kernels
	*/
	//@{
	Kernel<double> _d5s1K;		//!< Direction 5, scale 1
	Kernel<double> _d7s1K;		//!< Direction 7, scale 1
	Kernel<double> _d5s2K;		//!< Direction 5, scale 2
	Kernel<double> _d7s2K;		//!< Direction 7, scale 2
	//@}

#endif

	/** @name Cell activities
		Variables for all scales and orientations: most layers have 8 orientations and 2 scales.
	*/
	//@{
	double* _d0s1;		//!< Direction 0, scale 1
	double* _d1s1;		//!< Direction 1, scale 1
	double* _d2s1;		//!< Direction 2, scale 1
	double* _d3s1;		//!< Direction 3, scale 1
	double* _d4s1;		//!< Direction 4, scale 1
	double* _d5s1;		//!< Direction 5, scale 1
	double* _d6s1;		//!< Direction 6, scale 1
	double* _d7s1;		//!< Direction 7, scale 1

	double* _d0s2;		//!< Direction 0, scale 2
	double* _d1s2;		//!< Direction 1, scale 2
	double* _d2s2;		//!< Direction 2, scale 2
	double* _d3s2;		//!< Direction 3, scale 2
	double* _d4s2;		//!< Direction 4, scale 2
	double* _d5s2;		//!< Direction 5, scale 2
	double* _d6s2;		//!< Direction 6, scale 2
	double* _d7s2;		//!< Direction 7, scale 2
	//@}

	/** @name Storage files
	*/
	//@{
	std::ofstream _d0s1F;			//!< Stores data of m_d0s1;
	std::ofstream _d1s1F;			//!< Stores data of m_d1s1;
	std::ofstream _d2s1F;			//!< Stores data of m_d2s1;
	std::ofstream _d3s1F;			//!< Stores data of m_d3s1;
	std::ofstream _d4s1F;			//!< Stores data of m_d4s1;
	std::ofstream _d5s1F;			//!< Stores data of m_d5s1;
	std::ofstream _d6s1F;			//!< Stores data of m_d6s1;
	std::ofstream _d7s1F;			//!< Stores data of m_d7s1;

	std::ofstream _d0s2F;			//!< Stores data of m_d0s2;
	std::ofstream _d1s2F;			//!< Stores data of m_d1s2;
	std::ofstream _d2s2F;			//!< Stores data of m_d2s2;
	std::ofstream _d3s2F;			//!< Stores data of m_d3s2;
	std::ofstream _d4s2F;			//!< Stores data of m_d4s2;
	std::ofstream _d5s2F;			//!< Stores data of m_d5s2;
	std::ofstream _d6s2F;			//!< Stores data of m_d6s2;
	std::ofstream _d7s2F;			//!< Stores data of m_d7s2;
	//@}

	double _clip_thres;		//!< Constant to clip kernels
};
}
