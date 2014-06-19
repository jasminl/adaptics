/**
 * \file Level6.h
 * \brief 3D formotion motion segmentation
 */
#pragma once

#include "KernelLayer.h"
#include "Level5B.h"
#include "Input.h"

namespace visio
{

class CLevel5B;

/**
 * \brief Motion segmentation (3D formotion level6)
 */
class CLevel6: public KernelLayer
{

public:

	/**
	 *	Constructor
	*/
	CLevel6(int iter_per_frame, double dt,
			double A9, double D9, double C9,double B9,
			double A, double sigma,
			double sigmaP,
			unsigned int verticalSize, unsigned int horizontalSize,
			unsigned int totFrames,
			unsigned int nbDirections,
			double V, double sigmaV,
			double VTilde, double sigmaVTilde,
			CLevel5B* pL5,
			Input<int>* pInput, const std::string& output_folder,
			std::string name = "motionLevel6",
			bool record = false, Buffer* buffer = nullptr, bool testKernels = false);

	/**
	 *	Destructor
	 */
	~CLevel6(){};

	/**
	 *	Main computing method
	 */
	void compute(void* object);

	/**
	 *	Adjusts kernels so that they are removed if below threshold and modulated otherwise.
	 *	\param k Kernel to remove.
	 *	\param factor to use to modulate the kernel.
	 *	\param threshold Value under which to cut weights
	 */
	void prune_kernel(Kernel<double>& k, double factor, double threshold);

	/**
	 *	Attentional top down attention
	 *	\param x X coordinate of unit to modulate
	 *	\param y Y coordinate of unit to modulate
	 *	\param a Gaussian magnitude
	 *	\param sigma Spread of attentional Gaussian
	 *	\param attX X coordinate of center of Gaussian attentional top-down
	 *	\param attY Y coordinate of center of Gaussian attentional top-down
	 *	\return Strenght of attentional top-down
	 */
	double top_down_weight(unsigned int x, unsigned int y, double a, double sigma,
			double attX, double attY);

	/** @name Updating methods
	*/
	//@{
	/**
	 *	Update including convolution with kernel O but no depth inhibition
	 *	\param array Array to update
	 *	\param direction Direction (0-7) to update
	 *	\param pLayer Pointer to layer to update (for use with openmp)
	 */
	double* update_full_no_depth(double* array, unsigned int direction, CLevel6* pLayer);

	/**
	 *	Update including depth separation and Gaussian interactions
	 *	\param array Array to update
	 *	\param direction Direction (0-7) to update
	 *	\param pLayer Pointer to layer to update (for use with openmp)
	 */
	double* update_noO_depth_gaussian(double* array, unsigned int direction, CLevel6* pLayer);

	/**
	 *	Update including neither kernel O convolution nor depth inhibition, with Gaussian interactions
	 *	\param array Array to update
	 *	\param direction Direction (0-7) to update
	 *	\param pLayer Pointer to layer to update (for use with openmp)
	 */
	double* update_noO_noDepth_gaussian(double* array, unsigned int direction, CLevel6* pLayer);
	//@}

	/** @name previous updating methods
	 *	These functions would be used without Gaussian kernel depth suppression: currently they are not used
	 */
	//@{
	/**
	 *	Update including no kernel O convolution but with depth inhibition
	 */
	double* update_noO_depth(double* array, unsigned int direction, double* L5BInput,
			double* op, CLevel6* pLayer);

	/**
     *	Update including neither kernel O convolution nor depth inhibition
	 */
	double* update_noO_noDepth(double* array, unsigned int direction, double* L5BInput, CLevel6* pLayer);
	//@}

	/**
	 *	Overloaded method to print kernels to file.
	 */
	void print_kernels(const std::string& output_folder, bool recreate = false, string optName = "");

	/** @name Buffer management methods
	*/
	//@{
	void setBaseSB();
	void setBasePB();
	void setBaseV();
	void setBaseRL5();
	void setBaseRL6();
	//@}

private:

	/** @name Pre-calculated constants
	*/
	//@{
	double _cst27;		//!< 1-dt*A9
	double _cst28;		//!< -dt*A9*C9
	double _cst29;		//!< dt*A9;
	double _cst30;		//!< B9
	//@}

	std::vector<double> _dg;				//!< Directional gradient (w^De) same as CLayer5B
	Kernel<double> _P;						//!< Additional kernel for P
	Input<int>*	   _input;					//!< This provides the location of attentional focus

	std::vector<unsigned int>   _att_col;	//!< Column (x coord) of locus of attention
	std::vector<unsigned int>   _att_row;	//!< Row (y coord) of locus of attention

	std::vector<double> _V;					//!< Gaussian coefficients for cross-direction pooling of L5B signals
	std::vector<double> _V_tilde;			//!< Gaussian coefficients for cross-direction pooling of L6 scale 1 signals
};
}
