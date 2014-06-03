#pragma once
/*! \file parameters.h
    List of all parameters used in the simulations.
*/

/**	@name Time parameters 
	These parameters determine the temporal charaacteristics of the simulations.
*/
//@{
extern double g_dt;						//!< Euler integration step
extern double g_cyclePerSecond;			//!< Number of stimulus cycles per second
//@}

extern double g_ClipThreshold;			//!< Constant used to round values to zero (for numerical stability) 

/** @name Input dimensions
	These parameters determine the size of the input array, and their values should be defined in the input XML file.
*/
//@{
extern unsigned int g_verticalSize;
extern unsigned int g_horizontalSize;
//@}

extern unsigned int g_nbDirections;		//!< Number of motion directions (typically 8)

/** @name Kernel printing flags
	These parameters indicate whether to print a given level's kernels to a text file.
*/
//@{
extern bool g_printL3Kernels;
extern bool g_printL4Kernels;
extern bool g_printL5AKernels;
extern bool g_printL5BKernels;
extern bool g_printL6Kernels;
//@}

extern int  g_UpperMostLevel;			//!< Index which signifies up to which level to simulate

/** @name Recording flags
	These parameters indicate whether to record the activity of a given level to a (binary) file.
*/
//@{
extern bool   g_recordMotionL2;
extern bool   g_recordMotionL3;
extern bool   g_recordMotionL4;
extern bool   g_recordMotionL5A;
extern bool   g_recordMotionL5B;
extern bool   g_recordMotionL6;
//@}

/** @name Level 2 parameters
	These parameters are for the input stage (Level 2).
*/
//@{
extern double g_A1;
extern double g_B1;
extern double g_C1;
extern double g_A2;
extern double g_K2;

extern double g_A3;
extern double g_B3;
extern double g_C3;
extern double g_K3;

extern double g_A4;
extern double g_B4;
extern double g_C4;
extern double g_K4;
//@}

/** @name Level 3 parameters
	These parameters are for the short-range filter (Level 3).
*/
//@{
extern double g_A5;
extern double g_L3sigmaXscale1;
extern double g_L3sigmaYscale1;
extern double g_L3sigmaXscale2;
extern double g_L3sigmaYscale2;
extern double g_L3G;				
extern double g_L3theta1;
extern double g_L3theta2;
//@}

/** @name Level 4 parameters
	These parameters are for the spatial competition and opponent direction inhibition stage (Level 4).
*/
//@{
extern double g_A6;
extern double g_C6;
extern double g_D6;
extern double g_L4sigma1J;
extern double g_L4sigma2J;
extern double g_L4sigma1K;
extern double g_L4J;
extern double g_L4K;
//@}

/** @name Level 5A parameters
	These parameters are for the formotion selection stage (Level 5A).
*/
//@{
extern double g_A7;
extern double g_Ke;
extern double g_Kb;
extern double g_Kz;
extern double g_L5Asigma;
extern double g_L5I;
//@}

/** @name Level 5B parameters
	These parameters are for the long-range filter (Level 5B).
*/
//@{
extern double g_A8;
extern double g_D8;
extern double g_alpha;
extern double g_L5Bsigmaxs1;
extern double g_L5Bsigmays1;
extern double g_L5BsigmaPs1;
extern double g_L5Bsigmaxs2;
extern double g_L5Bsigmays2;
extern double g_L5BsigmaPs2;
extern double g_L5Ls1;
extern double g_L5thetas1;
extern double g_L5Ls2;
extern double g_L5thetas2;
extern double g_minW;
extern double g_maxW;
extern unsigned int g_wdeType;
//@}

/** @name Level 6 parameters
	These parameters are for the motion grouping stage (Level 6).
*/
//@{
extern double g_A9;
extern double g_D9;
extern double g_C9;
extern double g_L6sigmaP;
extern double g_V;
extern double g_sigmaV;
extern double g_VTilde;
extern double g_sigmaVTilde;
extern double g_B9;
//@}
