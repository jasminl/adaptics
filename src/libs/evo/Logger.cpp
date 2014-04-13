#include "evo/Logger.h"
#include <fstream>

using namespace std;

Logger::Logger(const string& folder, const string& file_name)
: _output_file(folder + "/" + file_name)
{}

void Logger::dump_stats(const tuple<int, float, float>& stats) const
{
	ofstream out(_output_file, ios_base::app);
	if(out.bad())
		throw ios_base::failure("Logger::dump_stats: unable to open " + _output_file + "\n");

	out<<get<0>(stats)<<" "<<get<1>(stats)<<" "<<get<2>(stats)<<endl;
}
