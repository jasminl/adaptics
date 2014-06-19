#include "evo/GrayCode.h"
#include "evo/Seq.h"
#include <boost/program_options.hpp>

using namespace std;
using namespace evo;

/**
 * Demonstrates the use of gray coding to decode sequences.
 */
int main(int argc, char* argv[])
{
	//Parse command line
	namespace po = boost::program_options;
	po::options_description desc("Allowed options");
	desc.add_options()
			("help", "graycodeapp, demonstrates the use of gray coding to encode sequences.")
			("len, l", po::value<int>()->required(), "Code length (i.e. number of bits in a gene)")
			("out, o", po::value<vector<string>>()->required(), "output file")
			("seq, s", po::value<int>()->required(), "Decoded sequence length");

	po::variables_map vm;
	po::store(po::parse_command_line(argc, argv, desc), vm);
	po::notify(vm);

	if (vm.count("help"))
	{
	    cout<<desc<<endl;
	    return 1;
	}

	int code_length  = vm["len"].as<int>();
	auto output_file = vm["out"].as<vector<string>>().front();

	int table_length = (int) pow((float)2.0, code_length);	//The number of sequences is a power of 2 of their length

	cout<<"File "<<output_file<<" will contain a table of "<<table_length<<" lines and "<<code_length<<" columns"<<endl;

	//Create Gray code table and output to file
	int** sequence = new int*[table_length];
	for(int i = 0; i < table_length; i++)
		sequence[i] = new int[code_length];
	auto table = create_gray(sequence, table_length, code_length);
	print_table(output_file, sequence, table_length, code_length);

	//Create a sequence and decode the genes
	int seq_len = vm["seq"].as<int>();
	Seq one_seq(seq_len, code_length);
	one_seq.randomize();

	//Find target "11" for example
	string target_tag(2, 0);
	target_tag[0] = 1;
	target_tag[1] = 1;
	auto instance = Seq::read_seq(one_seq.data(), seq_len, code_length, target_tag);

	//Show corresponding decoded sequence
	for(auto i: instance)
	{
		for(char p: i)
			cout<<(int)p;
		cout<<": "<<match_gray(i, table);
		cout<<endl;
	}
	//Cleanup
	for(int i = 0; i < table_length; i++)
		 delete[] sequence[i];
	delete[] sequence;

}
