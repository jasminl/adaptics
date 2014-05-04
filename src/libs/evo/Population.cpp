#include "evo/Population.h"
#include <stdexcept>
#include <cmath>
#include <string.h>

using namespace std;

namespace evo
{

Population::Population(int pop_size, int nb_building_blocks, int nb_neutral_blocks,
	float mutation_rate, float cross_over_rate, int building_block_size, int rr_type)
: _pop_size(pop_size), _nb_building_blocks(nb_building_blocks), _nb_neut_blocks(nb_neutral_blocks),
  _mutation_rate(mutation_rate), _cross_over_prob(cross_over_rate), _building_block_size(building_block_size),
  _gen_length(nb_building_blocks * (nb_neutral_blocks + building_block_size))
{
	//2D array of population data 
	_pop.resize(_pop_size);
	for(auto& pop: _pop)
		pop.resize(_gen_length);

	//1D array of fitness data
	_pop_fitness.resize(_pop_size, 0);

	//2D space for offspring generation
	_future_pop = _pop;

	//allocation of memory for royal road function
	switch(rr_type)
	{
	case FIXED_RR: //creates RR2
		_rr = new CRR2(nb_building_blocks, nb_neutral_blocks, building_block_size);
		break;
	case FLOATING_RR: //creates FRR2
		_rr = new CFRR2(nb_building_blocks, _gen_length, building_block_size);
		break;
	default:
		throw runtime_error("\nPopulation::Population: invalid population type");
	}
}

Population::~Population()
{
	delete _rr;
}

void Population::rnd_alloc()
{
	//fills the table for the entire population
	for(int i = 0; i < _pop_size; i++)
		for(int j = 0; j < _gen_length; j++)
			_pop[i][j] = char(rand() % 2);

	int rep;
	while(different != (rep = compare_seq()))
		sgl_alloc(rep);
}

void Population::sgl_alloc(int individual)
{
	if((individual < 0) || (individual > _pop_size))
		throw runtime_error("Population::sgl_alloc: attempted to replace individual beyond population range\n");

	for(int i = 0; i < _gen_length; i++)
		_pop[individual][i] = char(rand() % 2);
}

int Population::compare_seq()
{
	static int iterations = 0;	//to know how many times the function was called

	if(iterations > max_iterations)
		throw runtime_error("Population::compare_seq: algorithm does not converge\n");

	for(int i = 0; i < _pop_size - 1; i++)
	{
		for(int j = i + 1; j < _pop_size; j++)
		{
			int ss_sq = 0;
			char* p_seq0 = _pop[i].data();	//for capture by auto-vectorizer
			char* p_seq1 = _pop[j].data();
			for(int k = 0; k < _gen_length; k++)
				ss_sq += (p_seq0[k] - p_seq1[k]) * (p_seq0[k] - p_seq1[k]);

			if(ss_sq == 0)
				return j; //Need to change individual j
		}
	}
	return different;	//means every individual is different
}

void Population::cross_over(int first, int second, int offspring, int cross)
{
	const int cross_point = rand() % _gen_length;	//single crossover location

	//For capture by auto-vectorizer
	auto p_fut = _future_pop[offspring].data();
	auto p_par = _pop[first].data();
	auto p_par1 = _pop[second].data();
	const int gen_length = _gen_length;

	if(offspring == (_pop_size - 1))
	{
		if(cross == NO_MIX_PARENTS)
		{
			for(int d = 0; d < gen_length; d++)
				p_fut[d] = p_par[d];
		}
		else
		{
			for(int d = 0; d < cross_point; d++)
				p_fut[d] = p_par[d];

			for(int d = cross_point; d < gen_length; d++)
				p_fut[d] = p_par1[d];
		}

		return;
	}

	if(cross == NO_MIX_PARENTS)
	{	//no crossover applied
		memcpy(_future_pop[offspring].data(), _pop[first].data(), sizeof(char) * gen_length);
		memcpy(_future_pop[offspring + 1].data(), _pop[second].data(), sizeof(char) * gen_length);
		return;
	}

	auto p_fut1 = _future_pop[offspring + 1].data();	//For capture by auto-vectorizer

	//sequence beginning taken from parent first
	for(int i = 0; i < cross_point; i++)
	{
		p_fut[i] = p_par[i];
		p_fut1[i] = p_par1[i];
	}

	//sequence end taken from parent second
	for(int i = cross_point; i < gen_length; i++)
	{
		p_fut[i] = p_par1[i];
		p_fut1[i] = p_par[i];
	}
}

int Population::select_parent()
{
	static int seltime = 0, sel2 = 0;	//to know how many times the function was selected
	++sel2;	//keeps track of the number of selection made

	float pt = (_pop_size-1) * (rand() / (float)(RAND_MAX));	//to have a number between 0 and 1

	long double tempsum = 0;	//stores cumulative offspring

	bool inside=true;
	int i = 0;

	while(inside)
	{
		for(i = 0; i < _pop_size; i++)
		{
			tempsum += _pop_fitness[i];

			if(tempsum >= pt)
			{
				inside = false;
				break;
			}
			else
				inside = true;
		}
	}
	++seltime;
	return i;	//returns the selected individual
}

int Population::compute_pop_fitness()
{
	for(int i = 0; i < _pop_size; i++)
		_pop_fitness[i] = float(_rr->cmp_blcks(_pop[i].data()));

	//determines if an optimal individual was found: stops the experiment
	for(int i = 0; i < _pop_size; i++)
	{
		if(optimum_31 == _pop_fitness[i])
			return optim_indiv;	//means we have found an optimal individual
	}

	m_stats.compute_stats(_pop_fitness);

	//sigma scaling
	for(int i = 0;i < _pop_size; i++)
	{
		_pop_fitness[i] = 1 + ((_pop_fitness[i] - m_stats._average)/(2 * m_stats._stdev));

		if(_pop_fitness[i] > sigma_clip)
			_pop_fitness[i] = sigma_clip;
		else if(_pop_fitness[i] < 0.0)
			_pop_fitness[i] = float(scaling_fitness);
	}

	return continue_iterating;
}

int Population::generation()
{
	for(int j = 0; j < _pop_size; j++)	//clears the next generation's space
		for(int k = 0; k < _gen_length; k++)
			_future_pop[j][k] = 0;

	if(optim_indiv == compute_pop_fitness())
		return optim_indiv;

	int ofspnb = 0;	//number of offspring produced

	while(ofspnb <= _pop_size - 1)
	{
		int first = select_parent();	//first parent
		int second = select_parent();	//second parent

		if(first == second)
			continue;

		cross_over(first, second, ofspnb, cross_over_prob());	//crossover operator

		mutate(ofspnb);	//mutation operator

		if(ofspnb != (_pop_size - 1))
			mutate(ofspnb + 1);

		ofspnb += 2;
	}

	//transfer from m_futpop to m_pop
	_pop = _future_pop;

	return continue_iterating;
}

int Population::cross_over_prob()
{
	return (rand() / float(RAND_MAX) <= _cross_over_prob?MIX_PARENTS:NO_MIX_PARENTS);
}

void Population::mutate(int individual)
{
	for(int i = 0; i < _gen_length; i++)
	{
		if(rand() / (float)RAND_MAX <= _mutation_rate)
			_future_pop[individual][i] = !(_future_pop[individual][i]);
	}
}

tuple<int, float, float> Population::stats() const
{
	return tuple<int, float, float>{m_stats._sum, m_stats._average, m_stats._stdev};
}

}
