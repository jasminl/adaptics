#include "evo/CRR2.h"
#include <string.h>

using namespace std;

CRR2::CRR2(int nb_building_blocks, int nb_neut_blocks, int building_block_size)
: RoyalRoad(nb_building_blocks, building_block_size), _neutral_blck(nb_neut_blocks)
{}

int CRR2::cmp_blcks(const char* individual)
{
	memset(&_fitness[0], 0, sizeof(int) * _size);

	int i = 0;
	while (i <= (_nb_blocks - 1) * (_neutral_blck + _bld_blck_size))
	{
		for (int j = 0; j < _bld_blck_size; j++)
		{
			if (individual[i + j]
					!= _rr_value[i / (_neutral_blck + _bld_blck_size)][j])
				break;

			else if (j == _bld_blck_size - 1)
				_fitness[i / (_neutral_blck + _bld_blck_size)] = 1;
		}
		i += _neutral_blck + _bld_blck_size;
	}

	return fndupblck(); //checks for upper level building blocks
}
