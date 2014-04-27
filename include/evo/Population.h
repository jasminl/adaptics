#pragma once

#include <tuple>
#include "RoyalRoad.h"
#include "CRR2.h"
#include "CFRR2.h"
#include "Stats.h"

namespace evo
{

/**
 * Constants for population
 */
namespace nsPOP
{
	static const float sigma_clip = 1.5;
	static const float scaling_fitness = 0.1;	//!<sigma scaling of Tanese (1989)
};

/**
 * Generates and evolves a population
 */
class Population
{
public:
	//typedefs, etc.
	enum{FIXED_RR, FLOATING_RR};				//!<The type of royal road function used
	enum{CONTINUE = -2, OPTIMAL_INDIV = -1};	//!<Returned by the generation function to determine whether to continue or not
	enum{MIX_PARENTS, NO_MIX_PARENTS};			//!<Whether or not to mix parent dna during cross-over
public:

	/**
	 * Constructor creates a two dimensional array of char to store the entire population. Throws an exception if
	 * a wrong royal road type is specified.
	 * \param pop_size Population size
	 * \param nb_building_blocks Number of royal road blocks in genotype
	 * \param nb_neutral_blocks The number of neutral blocks per sequence
	 * \param mutation_rate The mutation rate
	 * \param cross_over_rate The cross over rate
	 * \param building_block_size The size of each building block (size 16 only supported so far)
	 * \param rr_type Type of royal road function to run
	 */
	Population(int pop_size, int nb_building_blocks, int nb_neutral_blocks,
		float mutation_rate, float cross_over_rate, int building_block_size, int rr_type);

	/**
	 * Frees memory allocated for the population
	 */
	virtual
	~Population();

	/**
	 * Fills up the dna sequence with rnd numbers. The usual rand() fct is used and srand()
	 * is used to set the starting value according to the actual time and to the number of
	 * time this function is called.
	 */
	void rnd_alloc();

	/**
	 * Fills in a specific individual's dna, generates a single DNA sequence for the individual specified
	 * by individual. Throws a runtime_error exception if the index specified is negative or greater than the
	 * population size.
	 * \param individual Index of the individual sequence to allocate.
	 */
	void sgl_alloc(int individual);

	/**
	 * Searches for duplicated sequences.
	 * \return If a duplicated sequence is found, the index of the duplicating sequence, DIFFERENT otherwise. Throws
	 * a runtime_error exception if the function is run beyond a specified number of maximum iterations.
	 */
	int compare_seq();

	/**
	 * Operates a single point crossover
	 * \param first The index of the first parent
	 * \param second The index of the second parent
	 * \param offspring The index of one of the potential two offsprings
	 * \param cross Whether or not to mix parent sequences
	 */
	void cross_over(int first, int second, int offspring, int cross = MIX_PARENTS);

	/**
	 * Chooses best fit parent according to roulette wheel selection
	 * \return The index of the selected parent
	 */
	int select_parent();

	/**
	 * Calculates fitness for entire population
	 * \return OPTIMAL_INDIV if an individual of theoretically optimal fitness is found, otherwise CONTINUE.
	 */
	int compute_pop_fitness();

	/**
	 * Performs one generation of evolution, computing fitness and applying cross-over and mutation operators.
	 * \return OPTIMAL_INDIV if an individual of theoretically optimal fitness is found, otherwise CONTINUE.
	 */
	int generation();

	/**
	 * Determines if crossover occurs according to probability
	 * \return MIX_PARENTS or NO_MIX_PARENTS depending on whether cross-over is determined to be applied.
	 */
	int cross_over_prob();

	/**
	 * Mutates the individual
	 * \param individual The index of the individual sequence to mutate.
	 */
	void mutate(int individual);

	/**
	 * Retrieve a copy of the stats structure for the current generation.
	 * \return A tuple of the sum, average and standard deviation of the population fitness
	 */
	std::tuple<int, float, float> stats() const;

protected:
	int _pop_size;				//!<The population size
	int _nb_building_blocks;	//!<The number of building blocks
	int _nb_neut_blocks;		//!<The number of neutral blocks
	float _mutation_rate;		//!<Mutation probability
	float _cross_over_prob;		//!<Crossover probability
	int _building_block_size;	//!<The size of building blocks
	const int _gen_length;			//!<Genotype length

	Stats m_stats;				//!<Stores simulation statistics

	RoyalRoad* _rr;				//!<The royal road signature

	std::vector<std::vector<char>> _pop;			//!<The current population
	std::vector<std::vector<char>> _future_pop;		//!<The next generation (necessary for crossover)
	std::vector<float> _pop_fitness;				//!<Array of population fitness
};

}
