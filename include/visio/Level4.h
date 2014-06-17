#pragma once

#include "KernelLayer.h"

namespace visio
{
/**
 * \brief Spatial competition and opponent direction of motion (3D formotion level 4)
 */
class CLevel4: public KernelLayer
{
public:
	/**
		Constructor
	*/
	CLevel4(int iter_per_frame, double dt,
		    double A6, double C6, double D6,
		    unsigned int verticalSize, unsigned int horizontalSize,
			unsigned int totFrames,
			unsigned int nbDirections,
		    double sigma1J, double sigma2J,
		    double sigma1K,
			double J, double K, const std::string& output_folder,
			std::string name = "motionLevel4",
			bool record = false,
			Buffer* buffer = nullptr,
			bool testKernels = false);

	/**
     *	Main computing method
	 */
	void compute(void* object);

	/**
	 *	Updates according to Eq. A11 from Berzh. et al. (2007), once the kernels have been convolved.
	 *	\param array Activity array to update.
	 *	\param bf1 First buffer storing result of convolutions.
	 *	\param bf2 Second buffer storing result of convolutions.
	 *	\param op Opposite direction activity.
	 *	\param direction Motion direction (0-7) being updated.
	 *	\return Pointer to updated activity.
	 */
	double* update(double* array, double* bf1, double* bf2, double* op, unsigned int direction);

	/**
	 *	Print kernels
	 */
	void print_kernels(const std::string& output_folder, bool recreate = false, std::string optName = "");

	/** Buffer management methods
	*/
	//@{
	void setBaseB2();
	void setBaseB();
	void setBaseSB();
	//@}

private:

	/**
	 * @name Pre-calculated constants
	 */
	//@{
	double _cst14;		//!< 1-dt*A6
	double _cst15;		//!< -dt*A6
	double _cst16;		//!< -dt*a6*C6
	double _cst17;		//!< dt*A6
	double _cst18;		//!< -dt*A6*0.1*C6
	double _cst19;		//!< -dt*A6*0.1*D6
	double _cst20;		//!< -dt*A6*D6
	//@}
};
}
