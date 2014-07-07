#include "evo/Population.h"
#include <boost/program_options.hpp>

#include <tuple>
#include <log4cxx/log4cxx.h>
#include <log4cxx/logger.h>
#include <log4cxx/basicconfigurator.h>
#include <log4cxx/fileappender.h>
#include <log4cxx/simplelayout.h>

using namespace std;
using namespace evo;

log4cxx::LoggerPtr logger(log4cxx::Logger::getLogger("sim"));	//Main logger

/**
 * Expects one argument that specify the log file.
 */
int main(int argc, char *argv[])
{
	//Parse command line
	namespace po = boost::program_options;
	po::options_description desc("Allowed options");
	desc.add_options()
			("help", "royalroadapp, demonstrates Royal Road functions for evolution.")
			("out, o", po::value<vector<string>>()->required(), "output file");

	po::variables_map vm;
	po::store(po::parse_command_line(argc, argv, desc), vm);
	po::notify(vm);

	if (vm.count("help"))
	{
		cout<<desc<<endl;
		return 1;
	}
	auto output_file = vm["out"].as<vector<string>>().front();

	//Setup logging to extneral file
	log4cxx::FileAppender* appender = new log4cxx::FileAppender(log4cxx::LayoutPtr(new log4cxx::SimpleLayout()),
			output_file, false);
	log4cxx::helpers::Pool p;
	appender->activateOptions(p);
	log4cxx::BasicConfigurator::configure(log4cxx::AppenderPtr(appender));

	//Init simulation
	int population_size = 2 * 2 * 256;
	int rr_fct = Population::FLOATING_RR;	//The type of royal road function
	int nb_building_blocks = 16;
	int nb_generations = 1000;
	int nb_neutral = 4;
	float mutation_rate = 0.005;
	float crossover_rate = 0.7;
	int bldblcksize = 2 * 8;

	Population pop(population_size, nb_building_blocks, nb_neutral,
			mutation_rate, crossover_rate, bldblcksize, rr_fct);

	//random allocation of each individual's DNA
	pop.rnd_alloc();

	for (int i = 0; i < nb_generations; i++)
	{
		int res = pop.generation();

		auto stat_info = pop.stats();
		LOG4CXX_INFO(logger, get<0>(stat_info)<<" "<<get<1>(stat_info)<<" "
				<<get<2>(stat_info)<<endl);

		if(Population::optim_indiv == res)
			break;	//Break if an optimal individual has been found
	}
	return 0;
}
