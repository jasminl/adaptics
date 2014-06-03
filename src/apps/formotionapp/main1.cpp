#include <iostream>
#include <iomanip>

#include <boost/program_options.hpp>
#include "visio/Layer.h"
#include "visio/Level2.h"
#include "visio/Level3.h"
#include "visio/Level4.h"
#include "visio/Level5A.h"
#include "visio/Level5B.h"
#include "visio/Level6.h"
#include "visio/Parameters.h"
#include "xmlUtils.h"

#include "DirManager.h"

#include <log4cxx/log4cxx.h>
#include <log4cxx/logger.h>
#include <log4cxx/basicconfigurator.h>

log4cxx::LoggerPtr logger(log4cxx::Logger::getLogger("sim"));	//Main logger

using namespace visio;
using namespace std;
int		     g_countDownToSave;			    //Indicates whether to save current iteration or not
 
int main(int argc, char *argv[], char *envp[])
{
	log4cxx::BasicConfigurator::configure();
	logger->setLevel(log4cxx::Level::getInfo());

	//Parse command line
	namespace po = boost::program_options;
	po::options_description desc("Allowed options");
	int initial_folder;
	desc.add_options()
			("help", "formotionapp, implements the 3D formotion simulations")
			("par, p", po::value<vector<string>>(), "parameter file/list")
			("in, i", po::value<vector<string>>(), "input file")
			("out, o", po::value<vector<string>>(), "output folder")
			("fol, f", po::value<int>(&initial_folder)->default_value(CDirManager::DEFAULT_START), "initial folder number")
			("sparse, s", "Whether to save all steps or the last one only")
			("l5", "Read in Layer5A values rather than compute them");

	po::variables_map vm;
	po::store(po::parse_command_line(argc, argv, desc), vm);
	po::notify(vm);

	if (vm.count("help"))
	{
	    cout<<desc<<endl;
	    return 1;
	}

	if (vm.count("par") == 0)
		throw runtime_error("Missing parameter file/folder");
	if(vm.count("in") == 0)
		throw runtime_error("Missing input file");
	if(vm.count("out") == 0)
		throw runtime_error("Missing output folder");

	bool save_all = true;
	if(vm.count("sparse"))
	{
		cout<<"Sparse saving configuration"<<endl;
		save_all = false;
	}

	bool load_layer5A = false;
	if(vm.count("l5"))
		load_layer5A = true;

	auto parameter = vm["par"].as<vector<string>>().front();
	auto input_path = vm["in"].as<vector<string>>().front();
	auto output_folder = vm["out"].as<vector<string>>().front();

	bool test_kernels = false;	//Indicates whether or not to calculate convolution times at each layer

	CDirManager d(output_folder, initial_folder);

	//Retrieve parameter and input file lists
	auto param_file = d.load_file(parameter);
	auto input_file = d.load_file(input_path);

	//Indicate how many files are parsed
	LOG4CXX_INFO(logger, "Setting up "<<param_file.size()
			<<" parameter configurations on "<<input_file.size()<<" stimuli.\n");
	try
	{
		for(auto p: param_file)
		{
			load_parameters(p);

			d.mk_par_dir();
			d.moveFile(p, (d.root_dir() + "/" + d.currentDir()).c_str());

			for(auto i: input_file)
			{
				//Input layer (Level 1)
				Input<int> input_layer(Input<int>::USE_MOTION, Input<int>::USE_V2, Input<int>::NO_USE_V1, Input<int>::USE_ATTEN);
				input_layer.import(i);

				//Determine number of cycles per second for current stimulus
				g_cyclePerSecond = input_layer.cycle_per_second(); //todo: shouldn't this be removed?

				//Determine number of iterations per frame
				unsigned int g_NbIterPerFrame = static_cast<unsigned int>(ceil(1/ (g_cyclePerSecond * input_layer.nb_frames() * g_dt)));
				LOG4CXX_INFO(logger,"Cycles/second: "<<g_cyclePerSecond<<endl
					<<"Number of frames: "<<input_layer.nb_frames()<<endl
					<<"Dt: "<<g_dt<<endl
					<<"Number of iterations per frame: "<<g_NbIterPerFrame<<endl);

				d.mk_stim_dir(i);
				d.moveFile(i, d.stim_dir().c_str());

				//Determine size of layers
				LOG4CXX_INFO(logger, "Height: "<<(g_verticalSize   = input_layer.height())
								<<"\tWidth: "<<(g_horizontalSize = input_layer.width())<<endl);

				Buffer Buffer(g_verticalSize, g_horizontalSize);

				//Declare motion layers
				CLevel2 motionLevel2(g_NbIterPerFrame, g_dt, g_A1, g_B1, g_C1, g_A2,
						g_K2, g_A3, g_B3, g_C3,
						g_K3, g_A4, g_B4, g_C4,
						g_K4, g_verticalSize, g_horizontalSize, g_nbDirections,
						input_layer.nb_frames(), d.stim_dir(), "motionLevel2", g_recordMotionL2);

				CLevel3 motionLevel3(g_NbIterPerFrame, g_dt, g_A5, g_L3G, g_verticalSize,
						g_horizontalSize, input_layer.nb_frames(), g_nbDirections, g_L3sigmaXscale1,
						g_L3sigmaYscale1, g_L3sigmaXscale2, g_L3sigmaYscale2, g_L3theta1,
						g_L3theta2, motionLevel2, d.stim_dir(), "motionLevel3",
						g_recordMotionL3, &Buffer, test_kernels);

				CLevel4 motionLevel4(g_NbIterPerFrame, g_dt, g_A6, g_C6, g_D6,
						g_verticalSize, g_horizontalSize, input_layer.nb_frames(), g_nbDirections,
						g_L4sigma1J, g_L4sigma2J, g_L4sigma1K, g_L4J,
						g_L4K, d.stim_dir(), "motionLevel4", g_recordMotionL4,
						&Buffer, test_kernels);

				CLevel5A motionLevel5A(g_NbIterPerFrame, g_dt, g_A7, g_Ke, g_Kb,
						g_Kz, g_verticalSize, g_horizontalSize, input_layer.nb_frames(),
						g_nbDirections, g_L5Asigma, g_L5I, input_layer.nb_V2_planes(),
						&input_layer, d.stim_dir(), "motionLevel5A", g_recordMotionL5A,
						&Buffer, test_kernels);

				CLevel5B motionLevel5B(g_NbIterPerFrame, g_dt, g_A8, g_D8, g_verticalSize,
						g_horizontalSize, input_layer.nb_frames(), g_nbDirections, g_alpha,
						g_L5Bsigmaxs1, g_L5Bsigmays1, g_L5BsigmaPs1, g_L5Ls1,
						g_L5thetas1, g_L5Bsigmaxs2, g_L5Bsigmays2, g_L5BsigmaPs2,
						g_L5Ls2, g_L5thetas2, g_minW, g_maxW,
						g_wdeType, d.stim_dir(), "motionLevel5B", g_recordMotionL5B,
						&Buffer, test_kernels);

				CLevel6 motionLevel6(g_NbIterPerFrame, g_dt, g_A9, g_D9, g_C9,
						g_B9, 1, NOT_ASSIGNED, g_L6sigmaP,
						g_verticalSize, g_horizontalSize, input_layer.nb_frames(),
						g_nbDirections, g_V, g_sigmaV, g_VTilde, g_sigmaVTilde,
						&motionLevel5B, &input_layer, d.stim_dir(), "motionLevel6",
						g_recordMotionL6, &Buffer, test_kernels);

				//Optionally print the kernels
				if(g_printL3Kernels)
					motionLevel3.print_kernels(d.stim_dir(), Kernel<double>::RECREATE_KER);
				if(g_printL4Kernels)
					motionLevel4.print_kernels(d.stim_dir(), Kernel<double>::RECREATE_KER);
				if(g_printL5AKernels)
					motionLevel5A.print_kernels(d.stim_dir(), Kernel<double>::RECREATE_KER);
				if(g_printL5BKernels)
					motionLevel5B.print_kernels(d.stim_dir(), Kernel<double>::RECREATE_KER);
				if(g_printL6Kernels)
					motionLevel6.print_kernels(d.stim_dir(), Kernel<double>::RECREATE_KER);

				//2. Iterate
				if(load_layer5A)
				{
					unsigned int nbf,  nbit;
					motionLevel5A.setL5ExternalInput(d.wkdir(), nbf, nbit);

					for(unsigned int f = 0; f < nbf; f++)
					{
						LOG4CXX_INFO(logger, "Processing frame "<<f<<endl);

						//Test whether to save all integration steps or only one
						if (!save_all)
							g_countDownToSave = g_NbIterPerFrame - 1;
						else
							g_countDownToSave = 0;

						for(unsigned int t = 0; t < nbit; t++)
						{
							if(g_UpperMostLevel > 5)
							{
								motionLevel5A.get_L5A_frame();
								motionLevel5B.compute(&motionLevel5A);
							}
							if(g_UpperMostLevel > 6)
								motionLevel6.compute(&motionLevel5B);

							if(!save_all)
								g_countDownToSave--;
						}
					}
				}
				else
				{
					for(unsigned int f = 0;f < input_layer.nb_frames(); f++)
					{
						LOG4CXX_INFO(logger, "Processing frame "<<f<<endl);

						input_layer.next_frame();

						//Test whether to save all integration steps or only one
						if (!save_all)
							g_countDownToSave = g_NbIterPerFrame - 1;
						else
							g_countDownToSave = 0;

						for(unsigned int t = 0; t < g_NbIterPerFrame; t++)
						{
							//For each timestep per frame
							if(g_UpperMostLevel > 1)
								motionLevel2.compute(input_layer.on_motion(), input_layer.off_motion());

							if(g_UpperMostLevel>2)
								motionLevel3.compute(&motionLevel2);

							if(g_UpperMostLevel>3)
								motionLevel4.compute(&motionLevel3);

							if(g_UpperMostLevel>4)
								motionLevel5A.compute(&motionLevel4);

							if(g_UpperMostLevel>5)
								motionLevel5B.compute(&motionLevel5A);

							if(g_UpperMostLevel>6)
								motionLevel6.compute(&motionLevel5B);

							if (!save_all)
								g_countDownToSave--;
						}
					}
				}

				LOG4CXX_INFO(logger, atoi(d.currentDir().c_str()));	//Indicate that simulation is done in log file
			}
		}
	}
	catch(exception& e)
	{	//Log exception
		LOG4CXX_INFO(logger, e.what());
		throw;
	}
	return 0;
}
