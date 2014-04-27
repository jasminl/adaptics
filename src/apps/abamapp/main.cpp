#include <iostream>
#include <vector>
#include "learn/Abam.h"

using namespace std;
using namespace learn;

int main(int argc, char *argv[], char *envp[])
{
	int nb_iterations = 100;	//Number of training iterations
	double dt = 0.01;		//Integration step
	double lrate = 0.2;		//Learning rate
	double maxd = 0.0001;	//Will stop settling iterations if difference below this threshold
	int maxit = 100;		//Max number of settling iterations

	//Two anti-correlated patterns
	vector<double> stis = {-1.0, 1.0, -1.0, 1.0, -1.0, 1.0};
	vector<double> stis2 = {1.0, -1.0, 1.0, -1.0, 1.0, -1.0};
	auto out = stis2;
	auto out2 = stis;

	Abam sim(dt, lrate, maxd, maxit, (int) stis.size(), (int)out.size());
	sim.init_matrix(Abam::RANDOM);

	for(int i = 0; i < nb_iterations; i++)
	{
		sim.train(stis, out);
		sim.train(stis2, out2);
	}

	cout<<"Error on first and second pattern: "<<sim.test(stis, out)<<" "<<sim.test(stis2, out2);

	return 0;
}
