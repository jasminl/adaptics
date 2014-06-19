/**
 * \file Level5B.h
 * \brief 3D formotion long-range filter
 */
#pragma once

#include "KernelLayer.h"
#include "Level6.h"

namespace visio
{

class CLevel6;

/**
 * \brief Long-range filter (3D formotion level 5B)
 */
class CLevel5B: public KernelLayer
{
public:

	/**
	 *	Constructor
	 */
	CLevel5B(int iter_per_frame, double dt,
		double A8, double D8,
		unsigned int verticalSize, unsigned int horizontalSize,
		unsigned int totFrames,
		unsigned int nbDirections,
		double alpha,
		double sigmaxs1, double sigmays1, double sigmaPs1, double Ls1, double thetas1,
		double sigmaxs2, double sigmays2, double sigmaPs2, double Ls2, double thetas2,
		double minW, double maxW, unsigned int type, const std::string& output_folder,
		std::string name = "motionLevel5B", bool record = false, Buffer* buffer = nullptr,
		bool testKernels = false);

	/** @name Buffer management methods
	*/
	//@{
	void setBaseB();
	void setBaseSB();
	void setBasePB();
	void setBaseV();
	//@}

	/**
	 *	Main computing method
	 */
	void compute(void* object);


	/**
	 * @name Updating methods
	 */
	//@{
#if NUM_THREADS==1
	/**
		Update function for non-separable kernels
	*/
	double* update(double* array, unsigned int direction,double* L5AInput, Kernel<double>& k, double theta,double* L6Selection);

	/**
		Update function for separable kernel
	*/
	double* updateSep(double* array, unsigned int direction,double* L5AInput, Kernel<double>& k, double theta,double* L6Selection);
#else

	/**
	 *	Update function for non-separable kernels
	 *	\param array Array to update
	 *	\param direction Direction of the update (0-7)
	 *	\param L5AInput Layer 5A boundary input
	 *	\param k Low-pass kernel applied to boundary input
	 *	\param pLayer Pointer to Level5B layer
	 *	\param b Buffer to store FFT output
	 *	\param preb Buffer to store rectified and squared input
	 *	\param theta Angle for update
	 *	\param L6Selection Layer 6 feedback input
	 *	\param scale (-1 or +1) which depth plane to update
	 */
	double* update(double* array, unsigned int direction, double* L5AInput, Kernel<double>& k,
			CLevel5B* pLayer, double* b, double* preb, double theta,
			double* L6Selection, int scale = -1);

	/**
	 *	Update function for separable kernel
	 */
	double* updateSep(double* array, unsigned int direction, double* L5AInput, Kernel<double>& k,
			CLevel5B* pLayer, double* b, double* sb, double* preb,
			double theta, double* L6Selection, int scale = -1);
#endif
	//@}

	/**
	 *	Overloaded method to print kernels.
	 */
	void print_kernels(const std::string& output_folder, bool recreate = false, std::string optName = "");

	/**
	 *	Sets pointer to CLevel6 object for feedback.
	 */
	void set_L6(CLevel6* object);

	/**
	 *	Access to the directional gradients.
	 *	\return Vector of directional gradients.
	 */
	std::vector<double> dg() const;

private:

	/**
	 * @name Thresholds
	 */
	//@{
	double _thetas1;		//!< For thresholding in Eq. A17, scale 1
	double _thetas2;		//!< For thresholding in Eq. A17, scale 2
	//@}

	/** @name Pre-calculated constants
	*/
	//@{
	double _cst25;			//!< 1-dt*A8
	double _cst26;			//!< dt*A8
	double _cst30;			//!< alpha*dt*A8
	//@}

	CLevel6* _L6;			//!< Pointer to Layer 6 for feedback

	std::vector<double> _dg;	//!< Directional gradient (w^De)

	/** @name feedback kernels
	*/
	//@{
	Kernel<double> _p_s1;		//!< Scale 1
	Kernel<double> _p_s2;		//!< Scale 2
	//@}
};

inline
void CLevel5B::set_L6(CLevel6* object)
{
	_L6 = object;
}

inline
vector<double> CLevel5B::dg() const
{
	return _dg;
}
}
