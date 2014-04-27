#pragma once

#include <string>
#include <vector>

namespace evo
{

/**
* Base class for royal road functions
*/
class RoyalRoad
{
public:

	/**
	 * Allocates all sequence-related buffers and generates the royal road functions. Throws a runtime_error
	 * exception in case the number of building blocks is not even.
	 * \param nb_building_blocks The number of building blocks in a sequence
	 * \param building_block_size The building block size
	 */
	RoyalRoad(int nb_building_blocks, int building_block_size);
	

	/**
	 * Destructor. No memory management performed here.
	 */
	virtual
	~RoyalRoad();

	/**
	 * Generates basic building blocks. Generates a runtime_error exception if it cannot create
	 * non-redundant blocks.
	 */
	void generate_rr();

	/**
	 * Ensures no two building blocks are the same
	 * \return DIFFERENT2 or NOTCONVERGE if it can or cannot find non-redundant sequences, respectively.
	 */
	int compare_seq();

	/**
	 * Replaces one block in the sequence
	 * \param block The index of the block to replace.
	 */
	void sgl_blck(int block);

	/**
	 * Computes the building block fitness. Must be overriden by all child classes.
	 */
	virtual
	int cmp_blcks(const char* individual) = 0;

	/**
	 * Calculates fitness according to a standard RR2 function. Only 16 blocks are supported for now.
	 * Throws a runtime_exception if the number of blocks is not 16.

     1. Basic Building blocks:

	 Block	Pattern   mem[]
	 -------------------------
	 0		********  0
	 1		********  1
	 2		********  2
	 ...    ********  ...
	 15		********  15
	 -------------------------
	 N=16
	 -------------------------

	 2. Upper level building blocks

	 Block	Component								mem[]
	 ------------------------------------------------------
	 16		0 1											16
	 17		2 3											17
	 18		4 5											18
	 19		6 7											19
	 20		8 9											20
	 21		10 11										21
	 22		12 13										22
	 23		14 15										23

	 24		0 1 2 3										24
	 25		4 5 6 7										25
	 26		8 9 10 11									26
	 27		12 13 14 15									27

	 28		0 1 2 3 4 5 6 7								28
	 29		8 9 10 11 12 13 14 15						29

	 30		0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15		30
	 ------------------------------------------------------
	 */
	int fndupblck();

protected:
	int _nb_blocks;								//!<Stores the number of basic building blocks
	int _bld_blck_size;							//!<Stores size of basic building blocks
	std::vector<int> _fitness;					//!<Array stores actual fitness
	int _size;									//!<Number of units in m_fit
	std::vector<std::vector<bool>> _rr_value;	//!<Contains the royal road function values
};

}
