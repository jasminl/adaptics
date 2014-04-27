#pragma once

#include <vector>

namespace evo
{

/**
 * Stores statistics for a given generation
 */
struct Stats
{
	Stats();

	/**
	 * Updates the stats
	 * \param fitness The entire population fitness vector
	 */
	void compute_stats(const std::vector<float>& fitness);

	int _sum;		//!<Total sum of all individual fitnesses
	float _average;	//!<Average fitness
	float _stdev;	//!<Standard deviation of fitness
};

}
