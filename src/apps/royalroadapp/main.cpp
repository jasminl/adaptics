#include "evo/Population.h"
#include "evo/constants.h"
#include "evo/Logger.h"

using namespace std;

int main(int argc, char *argv[])
{
	int population_size = 2 * 2 * 256;
	int rr_fct = Population::FLOATING_RR;	//The type of royal road function
	int nb_building_blocks = 16;
	int nb_generations = 1000;
	int nb_neutral = 4;
	float mutation_rate = 0.005;
	float crossover_rate = 0.7;
	int bldblcksize = 2 * 8;

	Logger log("/home/jasminl/workspace/cppmix/", "data");

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
