#include "evo/Stats.h"
#include <cmath>

using namespace std;

Stats::Stats()
: _sum(0), _average(0.0), _stdev(0.0)
{}

void Stats::compute_stats(const vector<float>& fitness)
{
	_sum = _stdev = _average = 0.0;

	for(int i = 0; i < (int) fitness.size(); i++)
		_sum += int(fitness[i]);

	_average = _sum / float(fitness.size());

	double ssquare = 0;
	for(int i = 0; i < (int)fitness.size(); i++)
		ssquare += pow((fitness[i] - _average), 2);

	_stdev = float(sqrt(ssquare/(fitness.size() - 1)));
}
