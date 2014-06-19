#include "evo/GrayCode.h"
#include <fstream>
#include <iostream>
#include <stdexcept>

using namespace std;

namespace evo
{

#define UNINIT_BIT -1	//!<Initial value to sequence bits to indicate they are not assigned values yet

int find_cur_pos(int* line, int code_length)
{
	int i = 0;
	for (i = 0; i < code_length; i++)
		if (*(line + i) != UNINIT_BIT)
			return i - 1;

	if (i >= code_length - 1)
		return code_length - 1;

	throw runtime_error("find_cur_pos: invalid position");
}

void inplace_relocate(int* dest, int* src, int start, int code_length)
{
	for (int i = start; i < code_length; i++)
		*(dest + i) = *(src + i);
}

map<string, int> create_gray(int** codes, int table_length, int code_length)
{
	//initialize the table to invalid codes
	for (int i = 0; i < table_length; i++)
		for (int j = 0; j < code_length; j++)
			*(*(codes + i) + j) = UNINIT_BIT;

	int nballoc = 1;
	while (nballoc < table_length)
	{
		nballoc *= 2; //when the thing is allocated

		int pos = find_cur_pos(codes[0], code_length);

		int i = 0;
		for (i = 0; i < (nballoc / 2); i++)
		{
			*(*(codes + i) + pos) = 0;
			inplace_relocate( *(codes + nballoc - i - 1), *(codes + i), pos + 1,
					code_length);
		}

		for (; i < nballoc; i++)
			*(*(codes + i) + pos) = 1;
	}

	//Fill map
	map<string, int> out;
	for(int c = 0; c < table_length; c++)
	{
		string cc(code_length, 0);
		for(int s = 0; s < code_length; s++)
			cc[s] = codes[c][s];
		out[cc] = c;
	}

	return out;
}

void print_table(const string& file_path, int** codes, int table_length, int code_length)
{
	ofstream file(file_path);
	if(file.bad())
		throw ios_base::failure("PrintTable: unable to open file for writing " + file_path);

	//Write header
	file<<"GrayCode";
	string buffer;
	buffer.resize(code_length - 9, ' ');
	file<<buffer<<"Value"<<endl;

	//Write sequences
	for (int i = 0; i < table_length; i++)
	{
		for (int j = 0; j < code_length; j++)
			file<<*(*(codes + i) + j);
		file<<"\t"<<i<<endl;
	}
}

int match_gray(const string& target, const map<string, int>& table)
{
	for(auto& pt: table)
	{
		if(pt.first.compare(target) == 0)
			return pt.second;
	}

	return gray_undefined;
}

}
