#pragma once

#include "KernelLayer.h"
#include "Level2.h"

namespace visio
{
/**
 * \brief Short-range filter (3D formotion level 3)
 */
class CLevel3: public KernelLayer
{
public:
	/**
	 *	Constructor
	 */
	CLevel3(int iter_per_frame, double dt, double A5, double G, unsigned int verticalSize, unsigned int horizontalSize,
			unsigned int totFrames, unsigned int nbDirections, double sigma1Scale1, double sigma2Scale1,
		    double sigma1Scale2, double sigma2Scale2, double theta1, double theta2,
			CLevel2& lev2, const std::string& output_folder, string name = "motionLevel3",
			bool record = false, Buffer* buffer = nullptr, bool testKernels = false);

	/**
	 *	Main compute method
	 */
	void compute(void* object);

	/** @name Buffer management methods
	*/
	//@{
	void setBaseB();
	void setBaseSB();
	void setBasePB();
	//@}

protected:	//Data members

	double _cst13;	//!< 1-A5*dt

	/** @name Thresholding parameters
	*/
	//@{
	double _theta1;
	double _theta2;
	//@}

	/** @name Indexing arrays
		Index arrays to speed up compute().
	*/
	//@{
	std::vector<Kernel<double>*> _ks1;
	std::vector<Kernel<double>*> _ks2;
	std::vector<std::string> _ns;
	std::vector<double*> _vs1;
	std::vector<double*> _vs2;
	//@}
};
}
