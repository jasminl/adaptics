#include "xmlUtils.h"
#include <iostream>
#include "tinyxml/tinyxml.h"
#include <fstream>
#include <stdexcept>
#include "visio/Parameters.h"

using namespace std;

void load_parameters(const string& filename)
{
	ifstream file(filename);
	if(file.bad())
		throw runtime_error("loadXMLParameters::file does not exist");

	TiXmlDocument doc(filename.c_str());
	if(!doc.LoadFile())
		throw ios_base::failure("loadXMLParameters::unable to parse file " + filename);
	else
	{
		TiXmlElement* root = doc.RootElement ();
		if(!root->NoChildren())
		{
			TiXmlNode* child = nullptr;
			while((child = root->IterateChildren(child)))
			{
				TiXmlElement* pChild = child->ToElement();//Convert to element to access attribute
				read_parameter(pChild->Attribute("name"), pChild->Attribute("value"));
			}
		}
		else
			throw runtime_error("loadXMLParameters: no parameter in file" + filename);
	}
}

void read_parameter(const char* vName, const char* vValue)
{
	if(!strcmp(vName,"g_dt"))
	{
		g_dt = atof(vValue);
		return;
	}

	if(!strcmp(vName,"g_nbDirections"))
	{
		g_nbDirections = atoi(vValue);
		return;
	}

	if(!strcmp(vName,"g_A1"))
	{
		g_A1 = atof(vValue);
		return;
	}

	if(!strcmp(vName,"g_B1"))
	{
		g_B1 = atof(vValue);
		return;
	}

	if(!strcmp(vName,"g_C1"))
	{
		g_C1 = atof(vValue);
		return;
	}

	if(!strcmp(vName,"g_A2"))
	{
		g_A2 = atof(vValue);
		return;
	}

	if(!strcmp(vName,"g_K2"))
	{
		g_K2 = atof(vValue);
		return;
	}

	if(!strcmp(vName,"g_A3"))
	{
		g_A3 = atof(vValue);
		return;
	}

	if(!strcmp(vName,"g_B3"))
	{
		g_B3 = atof(vValue);
		return;
	}

	if(!strcmp(vName,"g_C3"))
	{
		g_C3 = atof(vValue);
		return;
	}

	if(!strcmp(vName,"g_K3"))
	{
		g_K3 = atof(vValue);
		return;
	}

	if(!strcmp(vName,"g_A4"))
	{
		g_A4 = atof(vValue);
		return;
	}

	if(!strcmp(vName,"g_B4"))
	{
		g_B4 = atof(vValue);
		return;
	}

	if(!strcmp(vName,"g_C4"))
	{
		g_C4 = atof(vValue);
		return;
	}

	if(!strcmp(vName,"g_K4"))
	{
		g_K4 = atof(vValue);
		return;
	}

	if(!strcmp(vName,"g_ClipThreshold"))
	{
		g_ClipThreshold = atof(vValue);
		return;
	}

	if(!strcmp(vName,"g_A5"))
	{
		g_A5 = atof(vValue);
		return;
	}


	if(!strcmp(vName,"g_L3sigmaXscale1"))
	{
		g_L3sigmaXscale1 = atof(vValue);
		return;
	}

	if(!strcmp(vName,"g_L3sigmaYscale1"))
	{
		g_L3sigmaYscale1 = atof(vValue);
		return;
	}

	if(!strcmp(vName,"g_L3sigmaXscale2"))
	{
		g_L3sigmaXscale2 = atof(vValue);
		return;
	}

	if(!strcmp(vName,"g_L3sigmaYscale2"))
	{
		g_L3sigmaYscale2 = atof(vValue);
		return;
	}

	if(!strcmp(vName,"g_L3G"))
	{
		g_L3G = atof(vValue);
		return;
	}

	if(!strcmp(vName,"g_L3theta1"))
	{
		g_L3theta1 = atof(vValue);
		return;
	}

	if(!strcmp(vName,"g_L3theta2"))
	{
		g_L3theta2 = atof(vValue);
		return;
	}

	if(!strcmp(vName,"g_A6"))
	{
		g_A6 = atof(vValue);
		return;
	}


	if(!strcmp(vName,"g_C6"))
	{
		g_C6 = atof(vValue);
		return;
	}

	if(!strcmp(vName,"g_D6"))
	{
		g_D6 = atof(vValue);
		return;
	}

	if(!strcmp(vName,"g_L4sigma1J"))
	{
		g_L4sigma1J = atof(vValue);
		return;
	}

	if(!strcmp(vName,"g_L4sigma2J"))
	{
		g_L4sigma2J = atof(vValue);
		return;
	}

	if(!strcmp(vName,"g_L4sigma1K"))
	{
		g_L4sigma1K = atof(vValue);
		return;
	}

	if(!strcmp(vName,"g_L4J"))
	{
		g_L4J = atof(vValue);
		return;
	}

	if(!strcmp(vName,"g_L4K"))
	{
		g_L4K = atof(vValue);
		return;
	}


	if(!strcmp(vName,"g_A7"))
	{
		g_A7 = atof(vValue);
		return;
	}

	if(!strcmp(vName,"g_Ke"))
	{
		g_Ke = atof(vValue);
		return;
	}

	if(!strcmp(vName,"g_Kb"))
	{
		g_Kb = atof(vValue);
		return;
	}

	if(!strcmp(vName,"g_Kz"))
	{
		g_Kz = atof(vValue);
		return;
	}

	if(!strcmp(vName,"g_L5Asigma"))
	{
		g_L5Asigma = atof(vValue);
		return;
	}

	if(!strcmp(vName,"g_L5I"))
	{
		g_L5I = atof(vValue);
		return;
	}

	if(!strcmp(vName,"g_A8"))
	{
		g_A8 = atof(vValue);
		return;
	}

	if(!strcmp(vName,"g_D8"))
	{
		g_D8 = atof(vValue);
		return;
	}

	if(!strcmp(vName,"g_alpha"))
	{
		g_alpha = atof(vValue);
		return;
	}

	if(!strcmp(vName,"g_L5Bsigmaxs1"))
	{
		g_L5Bsigmaxs1 = atof(vValue);
		return;
	}

	if(!strcmp(vName,"g_L5Bsigmays1"))
	{
		g_L5Bsigmays1 = atof(vValue);
		return;
	}

	if(!strcmp(vName,"g_L5BsigmaPs1"))
	{
		g_L5BsigmaPs1 = atof(vValue);
		return;
	}

	if(!strcmp(vName,"g_L5Ls1"))
	{
		g_L5Ls1 = atof(vValue);
		return;
	}

	if(!strcmp(vName,"g_L5thetas1"))
	{
		g_L5thetas1 = atof(vValue);
		return;
	}

	if(!strcmp(vName,"g_L5Bsigmaxs2"))
	{
		g_L5Bsigmaxs2 = atof(vValue);
		return;
	}

	if(!strcmp(vName,"g_L5Bsigmays2"))
	{
		g_L5Bsigmays2 = atof(vValue);
		return;
	}

	if(!strcmp(vName,"g_L5BsigmaPs2"))
	{
		g_L5BsigmaPs2 = atof(vValue);
		return;
	}

	if(!strcmp(vName,"g_L5Ls2"))
	{
		g_L5Ls2 = atof(vValue);
		return;
	}

	if(!strcmp(vName,"g_L5thetas2"))
	{
		g_L5thetas2 = atof(vValue);
		return;
	}

	if(!strcmp(vName,"g_A9"))
	{
		g_A9 = atof(vValue);
		return;
	}

	if(!strcmp(vName,"g_D9"))
	{
		g_D9 = atof(vValue);
		return;
	}

	if(!strcmp(vName,"g_C9"))
	{
		g_C9 = atof(vValue);
		return;
	}
	if(!strcmp(vName,"g_L6sigmaP"))
	{
		g_L6sigmaP = atof(vValue);
		return;
	}

	if(!strcmp(vName,"g_minW"))
	{
		g_minW = atof(vValue);
		return;
	}

	if(!strcmp(vName,"g_maxW"))
	{
		g_maxW = atof(vValue);
		return;
	}

	if(!strcmp(vName,"g_wdeType"))
	{
		g_wdeType = atoi(vValue);
		return;
	}

	if(!strcmp(vName,"g_printL3Kernels"))
	{
		g_printL3Kernels = static_cast<bool>(atoi(vValue)!=0);
		return;
	}

	if(!strcmp(vName,"g_printL4Kernels"))
	{
		g_printL4Kernels = static_cast<bool>(atoi(vValue)!=0);
		return;
	}

	if(!strcmp(vName,"g_printL5AKernels"))
	{
		g_printL5AKernels = static_cast<bool>(atoi(vValue)!=0);
		return;
	}

	if(!strcmp(vName,"g_printL5BKernels"))
	{
		g_printL5BKernels = static_cast<bool>(atoi(vValue)!=0);
		return;
	}

	if(!strcmp(vName,"g_printL6Kernels"))
	{
		g_printL6Kernels = static_cast<bool>(atoi(vValue)!=0);
		return;
	}

	if(!strcmp(vName,"g_UpperMostLevel"))
	{
		g_UpperMostLevel = atoi(vValue);
		return;
	}

	if(!strcmp(vName,"g_recordMotionL2"))
	{
		g_recordMotionL2 = static_cast<bool>(atoi(vValue)!=0);
		return;
	}

	if(!strcmp(vName,"g_recordMotionL3"))
	{
		g_recordMotionL3 = static_cast<bool>(atoi(vValue)!=0);
		return;
	}

	if(!strcmp(vName,"g_recordMotionL4"))
	{
		g_recordMotionL4 = static_cast<bool>(atoi(vValue)!=0);
		return;
	}

	if(!strcmp(vName,"g_recordMotionL5A"))
	{
		g_recordMotionL5A = static_cast<bool>(atoi(vValue)!=0);
		return;
	}
	if(!strcmp(vName,"g_recordMotionL5B"))
	{
		g_recordMotionL5B = static_cast<bool>(atoi(vValue)!=0);
		return;
	}
	if(!strcmp(vName,"g_recordMotionL6"))
	{
		g_recordMotionL6 = static_cast<bool>(atoi(vValue)!=0);
		return;
	}
	if(!strcmp(vName,"g_V"))
	{
		g_V = atof(vValue);
		return;
	}
	if(!strcmp(vName,"g_sigmaV"))
	{
		g_sigmaV = atof(vValue);
		return;
	}
	if(!strcmp(vName,"g_VTilde"))
	{
		g_VTilde = atof(vValue);
		return;
	}
	if(!strcmp(vName,"g_sigmaVTilde"))
	{
		g_sigmaVTilde = atof(vValue);
		return;
	}
	if(!strcmp(vName,"g_B9"))
	{
		g_B9 = atof(vValue);
		return;
	}

	throw runtime_error("readParameters: unidentified parameter" + string(vName));
}

