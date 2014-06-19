/**
 * \file Buf.h
 * \brief Buffer class used by 3D formotion model
 *
 * Levels 3-6 for performing various operations such as convolutions and rectifications.
 * See 'bufferAssignments.xls' for a table detailing the use of each buffer.
 */
#pragma once

#include "config.h"

namespace visio
{

/**
 * \brief Buffers for openMP parallelism
 *
 * Contains a number of buffers which are accessed by model layers to perform various operations when using openMP.
 */
class Buffer
{
public:

	/**
		Constructor, allocates buffer space
		@param vsize the vertical size of displays processed.
		@param hsize the horizontal size of displays processed.
	*/
	Buffer(unsigned int vsize, unsigned int hsize);

	/**
		Destructor
	*/
	virtual
	~Buffer(void);

public:	//Data members

	unsigned int _vSize;	//!< Vertical size
	unsigned int _hSize;	//!< Horizontal size
	
	/**	@name Basic buffers
		These buffers are used for 1- 4- ans 16-threads versions.
	*/
	//@{
	double* _buffer0;
	double* _buffer1;
	double* _buffer2;
	double* _buffer3;
	double* _buffer4;
	double* _buffer5;
	double* _buffer6;
	double* _buffer7;
	double* _buffer8;
	double* _buffer9;
	double* _buffer10;
	double* _buffer11;
	double* _buffer12;
	double* _buffer13;
	double* _buffer14;
	double* _buffer15;
	double* _buffer16;
	double* _buffer17;
	double* _buffer18;
	double* _buffer19;
	double* _buffer20;
	double* _buffer21;
	double* _buffer22;
	double* _buffer23;
	double* _buffer24;
//@}

#if NUM_THREADS == 4 || NUM_THREADS == 16
	
	/**	@name Secondary buffers
		These buffers are used for 4- and 16-threads versions.
	*/
	//@{
	double* _buffer25;
	double* _buffer26;
	double* _buffer27;
	double* _buffer28;
	double* _buffer29;
	double* _buffer30;
	double* _buffer31;
	//@}
#endif

#if NUM_THREADS == 16

	/**	@name Tertiary buffers
		These buffers are used for 16-threads versions.
	*/
	//@{
	double* _buffer32;
	double* _buffer33;
	double* _buffer34;
	double* _buffer35;
	double* _buffer36;
	double* _buffer37;
	double* _buffer38;
	double* _buffer39;
	double* _buffer40;
	double* _buffer41;
	double* _buffer42;
	double* _buffer43;
	double* _buffer44;
	double* _buffer45;
	double* _buffer46;
	double* _buffer47;
	double* _buffer48;
	double* _buffer49;
	double* _buffer50;
	double* _buffer51;
	double* _buffer52;
	double* _buffer53;
	double* _buffer54;
	double* _buffer55;
	double* _buffer56;
	double* _buffer57;
	double* _buffer58;
	double* _buffer59;
	double* _buffer60;
	double* _buffer61;
	double* _buffer62;
	double* _buffer63;
	double* _buffer64;
	double* _buffer65;
	double* _buffer66;
	double* _buffer67;
	double* _buffer68;
	double* _buffer69;
	double* _buffer70;
	double* _buffer71;
	//@}
#endif

};
}
