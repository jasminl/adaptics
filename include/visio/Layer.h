#pragma once

#include <iostream>
#include <vector>
#include <map>
#include <list>

#include "visio/Buf.h"
#include "tinyxml/tinyxml.h"
#include "jama/jama_svd.h"
#include "tnt/tnt_array2d.h"
#include "tnt/tnt_array1d.h"
#include "fftw3.h"
#include "visio/config.h"

namespace visio
{
#define PI 3.14159265	//!< Value used for pi
#define THRES 0.000001	//!< Value below which certain quantities are considered to be zero
#define NOT_ASSIGNED -1	//!< Constant used to assign no particular scale to a given layer's kernel bank

/**
	Base class for all model levels.
*/	
class CLayer
{
public:
	/**
	 *	Constructor
	 *	\param VSize Vertical size of input frames.
	 *	\param HSize Horizontal size of input frames.
	 *	\param totFrames The total number of frames in the sequence.
	 *	\param nbDirections The number of motion directions to process (must be 8).
	 *	\param name The name of the layer.
	 *	\param record Indicates whether to record the activity of that level.
	 *	\param buffer Pointer to the buffer object use for convolutions and updating.
	 */
	CLayer(unsigned int VSize, unsigned int HSize, unsigned int totFrames, unsigned int nbDirections,
		  std::string name, bool record, Buffer* buffer);

	/**
	 *	Destructor
	 */
	virtual ~CLayer(void){}

	/**
	 *	Creates a 2D array.
	 *	\param pointer Pointer to array to be created.
	 *	\param vSize Vertical size of array.
	 *	\param hSize Horizontal size of array.
	 *	\return Pointer to the created array.
	 */
	double* create_array(double*& pointer, unsigned int vSize, unsigned int hSize);
	
	/**
	 *	Deletes a previously created array.
	 *	\param pointer Pointer to the array to delete.
	 */
	void delete_array(double*& pointer);

	/**
	 *	Shows activity of a given layer on stdout.
	 *	\param array Pointer to array to display.
	 */
	void show_activity(double* array) const;

	/**
	 *	Used to set the initial value of an array.
	 *	\param array Pointer to array to set the value of.
	 *	\param value Value at which to set each item in the array.
	 */
	void set_value(double* array, double value);

	/**
	 *	To access an array by name.
	 *	\param name Name of the array to access.
	 *	\return Pointer to the named array.
	 */
	virtual
	double* array(const char* name) = 0;

	/**
	 *	Indicates whether to record the activity of that layer or not.
	 *	\return Boolean flag indicating whether activity is recorded or not.
	 */
	bool keep_record() const;
 
	/**
	 *	Records the value of all arrays of that layer.
	 */
	virtual
	void save() = 0;

	/**
	 *	\return The vertical size of the frames.
	 */
	unsigned int vSize() const;

	/**
	 *	\return The horizontal size of the frames.
	 */
	unsigned int hSize() const;

	/**
	 *	Returns the offset that must be used to split in half an input layer (for subsequent parallelizing)
	 */
	unsigned int chop_input2(unsigned int vsize, unsigned int hsize);

	/**
	 *	Returns the offset that must be used to split in half an input layer (for subsequent parallelizing)
	 */
	static
	unsigned int chop_input2(unsigned int totSize);

	/** @name Buffer setting member functions
	 *	The following group of functions assigns IDs to various buffers, depending on the number of threads and the layer
 	 */
	//@{
	virtual void setBaseB(){};
	virtual void setBaseSB(){};
	virtual void setBasePB(){};
	virtual void setBaseB2(){};
	virtual void setBaseV(){};
	virtual void setBaseRL5(){};
	virtual void setBaseRL6(){};
	//@}

protected:

	unsigned int _vSize;			//!< Number of rows of units
	unsigned int _hSize;			//!< Number of columns
	unsigned int _nb_frames;		//!< The total number of frames in the simulation
	unsigned int _nb_dir;			//!< The number of directions

	std::string _level_name;		//!< Name of current level

	bool _record;					//!< Indicates whether or not to keep track of values in time

	Buffer* _buf;						//<! Unique buffer object, contains all buffers

	/** @name Buffers
	*/
	//@{
	std::vector<double*> _baseB;
	std::vector<double*> _baseSB;
	std::vector<double*> _basePB;
	std::vector<double*> _baseB2;
	std::vector<double*> _baseV;
	std::vector<double*> _baseRL5;
	std::vector<double*> _baseRL6;
	//@}
};

inline
bool CLayer::keep_record() const
{
	return _record;
}

inline
unsigned int CLayer::vSize() const
{
	return _vSize;
}

inline
unsigned int CLayer::hSize() const
{
	return _hSize;
}

}
