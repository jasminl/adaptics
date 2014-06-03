#pragma once

#include "Layer.h"
#include <fstream>

namespace visio
{
/**
	Level 2 (transient cells)
*/
class CLevel2: public CLayer
{
public:
	/**
     *	Base constructor
	 *	\sa CLayer
	 */
	CLevel2(int iter_per_frame, double dt,
		    double A1,double B1,double C1,double A2,double K2,
			double A3,double B3,double C3,double K3,
			double A4,double B4,double C4,double K4,
			unsigned int verticalSize, unsigned int horizontalSize,
			unsigned int nbDirections, unsigned int totFrames,
			const std::string& output_folder, string name = "motionLevel2",
			bool record = false);

	/**
	 *	Destructor
 	 */
	virtual ~CLevel2();

public: //Member functions

	/**
		Main calculations for the transient cell level.
		\param motionOn Input to the ON motion stream.
		\param motionOff Input to the OFF motion stream.
	*/
	void compute(const std::vector<int>& motionOn, const std::vector<int>& motionOff);

	/**
		Returns an array by name.
		\param name Name of the array to return.
		\return A pointer to the desired array.
	*/
	double* array(const char* name);

	/** @name Access methods.
		Methods to access various arrays directly.
	*/
	//@{
	double* xon(){return _xOn;}
	double* xoff(){return _xOff;}
	double* zxon(){return _zOn;}
	double* zoff(){return _zOff;}
	double* b(){return _b;}
	double* c0(){return _c0;}
	double* c1(){return _c1;}
	double* c2(){return _c2;}
	double* c3(){return _c3;}
	double* c4(){return _c4;}
	double* c5(){return _c5;}
	double* c6(){return _c6;}
	double* c7(){return _c7;}
	double* e0(){return _e0;}
	double* e1(){return _e1;}
	double* e2(){return _e2;}
	double* e3(){return _e3;}
	double* e4(){return _e4;}
	double* e5(){return _e5;}
	double* e6(){return _e6;}
	double* e7(){return _e7;}
	//@}

	/**
		Open output files to store results.
		\param output_folder Folder in which to open files
		\param iter_per_frame Number of numerical integration steps per stimulus frame
		\param dt Integration step
	*/
	void open_files(const std::string& output_folder, int iter_per_frame, double dt);

	/**
		Save data to files.
	*/
	void save();

protected: //Data members

	/**
	 * @name Pre-calculated constants
	 */
	//@{
	//Constants for X (Eq.A5)
	double _cst1;	//!< 1- dt*A1*B1
	double _cst2;  //!< dt*A1*C1
	double _cst3;  //!< -dt*A1

	//Constants for Z (Eq.A6)
	double _cst4;  //!< dt*A2
	double _cst5;  //!< 1-dt*A2
	double _cst6;	//!< -dt*A2*K2

	//Constants for C (Eq.A7)
	double _cst7;	//!< 1-dt*A3*B3
	double _cst8;	//!< dt*A3*C3
	double _cst9;  //!< -dt*A3*K3

	//Constants for E (Eq.A8)
	double _cst10;	//!< 1 - dt*A4*B4
	double _cst11; //!< dt*A4*C4
	double _cst12; //!< -dt*A4*K4
	//@}

	/** @name Variable arrays
	*/
	//@{
	double* _xOn;
	double* _xOff;

	double* _zOn;
	double* _zOff;

	double* _b;

	double* _c0;
	double* _c1;
	double* _c2;
	double* _c3;
	double* _c4;
	double* _c5;
	double* _c6;
	double* _c7;

	double* _e0;
	double* _e1;
	double* _e2;
	double* _e3;
	double* _e4;
	double* _e5;
	double* _e6;
	double* _e7;
	//@}

	/** @name Storing variable arrays
	*/
	//@{
	std::ofstream  _xOnStoreF;
	std::ofstream  _xOffStoreF;

	std::ofstream  _zOnStoreF;
	std::ofstream  _zOffStoreF;

	std::ofstream  _bStoreF;

	std::ofstream  _c0StoreF;
	std::ofstream  _c1StoreF;
	std::ofstream  _c2StoreF;
	std::ofstream  _c3StoreF;
	std::ofstream  _c4StoreF;
	std::ofstream  _c5StoreF;
	std::ofstream  _c6StoreF;
	std::ofstream  _c7StoreF;

	std::ofstream  _e0StoreF;
	std::ofstream  _e1StoreF;
	std::ofstream  _e2StoreF;
	std::ofstream  _e3StoreF;
	std::ofstream  _e4StoreF;
	std::ofstream  _e5StoreF;
	std::ofstream  _e6StoreF;
	std::ofstream  _e7StoreF;
	//@}
};
}
