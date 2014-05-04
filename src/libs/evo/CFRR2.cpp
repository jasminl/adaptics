#include "evo/CFRR2.h"
#include <string.h>

using namespace std;

namespace evo
{

CFRR2::CFRR2(int NbBldBlck, int dnalength, int BldBlckSze)
: RoyalRoad(NbBldBlck, BldBlckSze), _dna_length(dnalength), _start_tag({true, true})
{}

int CFRR2::cmp_blcks(const char* individual)
{
	memset(&_fitness[0], 0, _size * sizeof(int));

	int strpos = 0, blckpos;

	while (strpos < _dna_length)
	{
		blckpos = find_tag(individual, strpos); //finds the start tag

		if (strpos == not_found)
			return not_found;

		find_block(individual, blckpos); //finds occurence of bld blcks
		++strpos;
	}

	return fndupblck();
}

int CFRR2::find_tag(const char* individual, int& position)
{
	for (int i = position; i < _dna_length - (_bld_blck_size + 1); i++)
	{
		if ((individual[i] == _start_tag[0]) && (individual[i + 1] == _start_tag[1]))
		{
			position = i;
			return i + 2; //starting position of the building blck
		}
	}

	return not_found; //means the start tag is not in the sequence from
					 //position to the end
}

void CFRR2::find_block(const char* individual, int position)
{
	for (int i = 0; i < _nb_blocks; i++)
	{
		int j;
		for (j = 0; j < _bld_blck_size; j++)
			if (_rr_value[i][j] != individual[position + j])
				break;

		if ((j == _bld_blck_size) && (_fitness[i] == 0))
			_fitness[i] = 1;
	}
}

}
