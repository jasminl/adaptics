#pragma once

#include <vector>
#include "RoyalRoad.h"

namespace evo
{

/**
 * \brief A floating representation version of CRR2
 */
class CFRR2: public RoyalRoad
{
public:
	//typedefs etc.
	static const int not_found = -3;	//!<Means no building blocks were found
public:

	/**
	 * Constructor, performs no other operations than calling the base class constructor and assigning
	 * the sequence length and start tag.
	 * \param nb_building_blocks The number of building blocks in a sequence
	 * \param dna_length Sequence length
	 * \param building_block_size The size of each building block
	 */
	CFRR2(int nb_building_blocks, int dna_length, int building_block_size);

	/**
	 * Finds instances of building blocks for floating representation and returns the total fitness. Throws
	 * a runtime_error exception in case of invalid block size.
	 * \param individual Const pointer to the individual sequence for which to seek royal road blocks
	 * \return The overall fitness of the sequence, or NOTFOUND if no start tag is present
	 */
	int cmp_blcks(const char* individual);

	/**
	 * Finds the start tag and returns a reference to its position
	 * \param individual Const pointer to the individual sequence for which to seek a start tag
	 * \param position Reference-passed argument into which the start tag's position is returned.
	 * \return The position of the start tag, or NOTFOUND if no start tag is found.
	 */
	int find_tag(const char* individual, int& position);

	/**
	 * Determines if a building block is present at the specified position. Increments the fitness by 1
	 * if the building block is found.
	 */
	void find_block(const char* individual, int position);

private:
	int _dna_length;				//!<The sequence length
	std::vector<bool> _start_tag;	//!<Contains start tag for building blocks
};

}
