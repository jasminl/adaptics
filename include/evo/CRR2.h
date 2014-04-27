#pragma once

#include "RoyalRoad.h"

namespace evo
{

/**
 * A fixed representation Royal Road functions
 */
class CRR2: public RoyalRoad
{
public:

	/**
	 * Constructor, performs no other operations than calling the base class constructor and assigning
	 * the number of neutral blocks.
	 * \param nb_building_blocks The number of building blocks in a sequence
	 * \param nb_neut_blocks The number of neutral blocks
	 * \param building_block_size The size of each building block
	 */
	CRR2(int nb_building_blocks, int nb_neut_blocks, int building_block_size);

	/**
	 * Finds instances of building blocks for fixed representation and returns the total fitness. Throws
	 * a runtime_error exception in case of invalid block size.
	 * \param individual Const pointer to the individual sequence for which to seek royal road blocks
	 * \return The overall fitness of the sequence
	 */
	int cmp_blcks(const char* individual);

private:
	int _neutral_blck;		//!<Number of neutral blocks (skipped by evaluation function)
};

}
