#include "evo/Population.h"
#include "evo/constants.h"
#include "evo/Logger.h"
#include <stdexcept>

using namespace std;
using namespace evo;

/**
 * Expects two arguments: a folder and a filename in which to store the logger's results
 */
int main(int argc, char *argv[])
{
	if(argc < 3)
		throw runtime_error("royalroadapp::main: must specify a folder and file in separate arguments");

	int population_size = 2 * 2 * 256;
	int rr_fct = Population::FLOATING_RR;	//The type of royal road function
	int nb_building_blocks = 16;
	int nb_generations = 1000;
	int nb_neutral = 4;
	float mutation_rate = 0.005;
	float crossover_rate = 0.7;
	int bldblcksize = 2 * 8;

	Logger log(argv[1], argv[2]);

	Population pop(population_size, nb_building_blocks, nb_neutral,
			mutation_rate, crossover_rate, bldblcksize, rr_fct);

	//random allocation of each individual's DNA
	pop.rnd_alloc();

	for (int i = 0; i < nb_generations; i++)
	{
		int res = pop.generation();
		log.dump_stats(pop.stats());

		if(OPTIMAL_INDIV == res)
			break;	//Break if an optimal individual has been found
	}
	return 0;
}
