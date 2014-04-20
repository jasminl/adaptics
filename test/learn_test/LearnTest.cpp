#define BOOST_TEST_MODULE Main
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include "learn/Abam.h"

using namespace std;

BOOST_AUTO_TEST_SUITE (test_suite_learn)

/**
 * Abam forward propagation
 */
BOOST_AUTO_TEST_CASE (AbamForward)
{
	double dt = 0.01;		//Integration step
	double lrate = 0.2;		//Learning rate
	double maxd = 0.0001;	//Will stop settling iterations if difference below this threshold
	int maxit = 100;		//Max number of settling iterations

	//Two anti-correlated patterns
	vector<double> stis = {1.0, 1.0};
	auto out = stis;

	Abam sim(dt, lrate, maxd, maxit, (int) stis.size(), (int)out.size());
	sim.init_matrix(Abam::ONES);
	sim.forward(stis);

	auto output = sim.output_layer();

	BOOST_CHECK_EQUAL(output.size(), 2);		//Check correct output layer size
	BOOST_CHECK_EQUAL(output[0], output[1]);	//Check same output at each site
	BOOST_CHECK_EQUAL(output[0], 1);
}

BOOST_AUTO_TEST_SUITE_END ()
