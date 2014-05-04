#include "evo/RoyalRoad.h"
#include <iostream>
#include <fstream>
#include <cmath>
#include <sstream>
#include <stdexcept>

using namespace std;

namespace evo
{

RoyalRoad::RoyalRoad(int nb_building_blocks, int building_block_size)
: _nb_blocks(nb_building_blocks), _bld_blck_size(building_block_size)
{
	if (nb_building_blocks % 2) //the number of basic building bloc must be even
		throw runtime_error("RoyalRoad::RoyalRoad: Uneven number of basic building blocks\n");

	/*
	 * Creates the table of fitness for a RR. The number of units in the table is equal
	 * to the total number of building blocks (including base and higher order)
	 */
	_size = 0;
	int i = 0;
	do
	{
		_size += int(nb_building_blocks / pow(2, i++));
	}
	while ((int(nb_building_blocks / pow(2, i - 1))) != 1);

	_fitness.resize(_size, 0);

	//now we create the royal road segments;
	_rr_value.resize(nb_building_blocks);
	for(auto& p_rr: _rr_value)
		p_rr.resize(building_block_size);

	generate_rr(); //generates the rr function
}

RoyalRoad::~RoyalRoad()
{}

void RoyalRoad::generate_rr()
{
	for (int i = 0; i < _nb_blocks; i++)
		sgl_blck(i);

	int rep = 0;
	while (rep != different)
	{
		rep = compare_seq();

		if (rep == different)
			break;

		if (not_converge == rep)
		{
			stringstream msg;
			msg<<"Cannot generate " << _nb_blocks	<< " different building blocks." << endl;
			throw runtime_error("RoyalRoad::generate_rr: " + msg.str() + "\n");
		}

		sgl_blck(rep); //to replace the faulty building block
	}
}

void RoyalRoad::sgl_blck(const int blck)
{
	for (int i = 0; i < _bld_blck_size; i++)
		_rr_value[blck][i] = bool(rand() % 2);
}

int RoyalRoad::compare_seq()
{
	static int iterations = 0;

	if (iterations > max_iterations)
		return not_converge; //means the algorithm does not converge

	for (int i = 0; i < _nb_blocks - 1; i++)
	{
		for (int j = i + 1; j < _nb_blocks; j++)
		{
			int k;
			for (k = 0; k < _bld_blck_size; k++)
			{
				if (_rr_value[i][k] != _rr_value[j][k])
					break;
			}
			if ((k == _bld_blck_size) && (_rr_value[i][k] == _rr_value[j][k]))
			{
				iterations++;
				return j; //means we want to change individual j
			}
		}
	}

	return different; //means every individual is different
}

int RoyalRoad::fndupblck()
{
	if (_nb_blocks != 16)
		throw runtime_error("RoyalRoad::fndupblck: invalid number of blocks\n");

	for (int i = 0, j = 0; i < 15; i += 2, j++)
		if (_fitness[i] && _fitness[i + 1])
			_fitness[16 + j] = 1;

	for (int i = 16, j = 0; i < 23; i += 2, j++)
		if (_fitness[i] && _fitness[i + 1])
			_fitness[24 + j] = 1;

	for (int i = 24, j = 0; i < 27; i += 2, j++)
		if (_fitness[i] && _fitness[i + 1])
			_fitness[28 + j] = 1;

	if (_fitness[28] && _fitness[29])
		_fitness[30] = 1;

	int partsum = 0;

	for (int i = 0; i < 31; i++)
		partsum += _fitness[i];

	return partsum;
}

}
