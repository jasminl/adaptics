#include <evo/Seq.h>

using namespace std;

namespace evo
{

Seq::Seq(int seq_len, int gene_len)
: _seq_len(seq_len), _gene_len(gene_len), _data(nullptr)
{}

Seq::~Seq()
{
	delete[] _data;
}

void Seq::inplace_int2string(int** codes, int table_length, int code_length,
		char** string_table)
{
	string_table = new char*[table_length];

	for (int i = 0; i < table_length; i++)
		*(string_table + i) = new char[code_length];

	for (int i = 0; i < table_length; i++)
	{
		for (int j = 0; j < code_length; j++)
			*(*(string_table + i) + j) = *(*(codes + i) + j);
	}
}

list<string> Seq::read_seq(const char* str, int str_len, int gene_length, const std::string& tag)
{
	int position = 0;
	string cppstr(str, str_len);

	list<string> out;
	while((int)string::npos != (position = cppstr.find(tag, position)))
	{
		if(position < str_len - gene_length - (int)tag.size())
			out.push_back(cppstr.substr(position + (int)tag.size(), gene_length));
		++position;	//Advance one position to ensure not decoding the same location
	}
	return out;
}

void Seq::randomize()
{
	if(_data != nullptr)
		delete[] _data;
	_data = new char[_seq_len];
	for(int i = 0; i < _seq_len; i++)
		_data[i] = rand() % 2;
}
}
