#include "visio/Parameters.h"

//This contains global variables

double g_dt;
double g_cyclePerSecond;
double g_ClipThreshold;

//Input dimensions
unsigned int g_verticalSize;
unsigned int g_horizontalSize;

//Number of directions
unsigned int g_nbDirections;

//Misc parameters
bool g_printL3Kernels;
bool g_printL4Kernels;
bool g_printL5AKernels;
bool g_printL5BKernels;
bool g_printL6Kernels;
int  g_UpperMostLevel;

bool g_recordMotionL2;
bool g_recordMotionL3;
bool g_recordMotionL4;
bool g_recordMotionL5A;
bool g_recordMotionL5B;
bool g_recordMotionL6;

//Level 2 parameters

double g_A1;
double g_B1;
double g_C1;
double g_A2;
double g_K2;

double g_A3;
double g_B3;
double g_C3;
double g_K3;

double g_A4;
double g_B4;
double g_C4;
double g_K4;

//Level 3 parameters
double g_A5;
double g_L3sigmaXscale1;
double g_L3sigmaYscale1;
double g_L3sigmaXscale2;
double g_L3sigmaYscale2;
double g_L3G;		
double g_L3theta1;
double g_L3theta2;

//Level 4 parameters
double g_A6;
double g_C6;
double g_D6;
double g_L4sigma1J;
double g_L4sigma2J;
double g_L4sigma1K;
double g_L4J;
double g_L4K;

//Level 5A parameters
double g_A7;
double g_Ke;
double g_Kb;
double g_Kz;
double g_L5Asigma;
double g_L5I;

//Level 5B parameters
double g_A8;
double g_D8;
double g_alpha;
double g_L5Bsigmaxs1;
double g_L5Bsigmays1;
double g_L5BsigmaPs1;
double g_L5Bsigmaxs2;
double g_L5Bsigmays2;
double g_L5BsigmaPs2;
double g_L5Ls1;
double g_L5thetas1;
double g_L5Ls2;
double g_L5thetas2;
double g_minW;
double g_maxW;
unsigned int g_wdeType;

//Level 6 parameters
double g_A9;
double g_D9;
double g_C9;
double g_L6sigmaP;
double g_V;
double g_sigmaV;
double g_VTilde;
double g_sigmaVTilde;
double g_B9;
