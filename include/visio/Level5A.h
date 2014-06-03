#pragma once

#include "KernelLayer.h"
#include "Input.h"

namespace visio
{
/**
	CLevel5A (formotion selection)
*/
class CLevel5A: public KernelLayer
{
public:
	/**
	 *	Constructor
	 */
	CLevel5A(int iter_per_frame, double dt,
			 double A7, double Ke, double Kb, double Kz,
			 unsigned int verticalSize, unsigned int horizontalSize,
			 unsigned int totFrames,
			 unsigned int nbDirections,
		     double sigma,
			 double I,
			 unsigned int nbV2Depth,
			 Input<int>* inputLayer, const std::string& output_folder,
			 std::string name = "motionLevel5A",
			 bool record = false,
			 Buffer* buffer = nullptr,
			 bool testKernels = false);

	/**
	 *	Destructor
 	 */
	~CLevel5A();

	/**
	 *	Main computing method
	 */
	void compute(void* object);

	/**
	 *	Performs Eq. A14 from Berzh. et al. 2007.
	 *	\param array Activity array to update.
	 *	\param L4Input Input from Level 4.
	 *	\param V2Input V2 boundary input.
	 *	\param bf1 Buffer to use for calculations.
	 *	\return A pointer to the array of updated activity.
	 */
	double* update(double* array, double* L4Input, double* V2Input, double* bf1);

	/**
	 *	Overloaded method to print kernels to file.
	 */
	void print_kernels(const std::string& output_folder, bool recreate = false, std::string optName = "");

	/**
	 *	Used if not simulating the front-end before the long-range filter but rather load stored data
	 *	from the previous layer.
	 */
	void setL5ExternalInput(const char* wd, unsigned int& nbFrames, unsigned int& nbIterations,
			std::vector<std::string> fname = std::vector<std::string>());

	/**
	 *	Obtain one frame from the formotion selection level to input to the long-range filter.
	 */
	void get_L5A_frame();

	/**
	 * @name Buffer management methods
	 */
	//@{
	void setBaseB();
	void setBaseSB();
	//@}

private:

	/** @name Pre-calculated constants
	*/
	//@{
	double _cst21;			//!< 1-dt*A7
	double _cst22;			//!< dt*A7*Ke
	double _cst23;			//!< dt*A7*Kz
	double _cst24;			//!< -dt*A7*Kb
	//@}

	Input<int>* _input;		//!< The input layer from which to get V2 boundaries

	double** _V2;			//!< Current V2 boundary input from _input
	unsigned int _nb_depth;	//!< The number of V2Depth planes

	std::vector<std::ifstream*> _input_file;	//!< Stored Layer 5A values
};
}
