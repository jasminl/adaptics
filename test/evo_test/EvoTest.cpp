#define BOOST_TEST_MODULE Main
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include "../src/libs/evo/ESDna.h"

using namespace std;
using namespace evo;

BOOST_AUTO_TEST_SUITE (test_suite_learn)

/**
 * Tests construction of ESDNa for type double
 */
BOOST_AUTO_TEST_CASE(ESctor_double)
{
	int len = 10;
	ESDna<double, double> dna(len);
	BOOST_CHECK_EQUAL(dna.len(), len);
}

/**
 * Test constructor with variable assigment
 */
BOOST_AUTO_TEST_CASE(ESctor_assign)
{
	int len = 10;
	vector<double> dna(len);
	for(int i = 0 ; i < len; i++)
		dna[i] = double(i);
	vector<double> mut = dna;

	ESDna<double, double> individual(dna, mut);

	for(int i = 0; i < len; i++)
	{
		BOOST_CHECK_EQUAL(dna[i], individual.dna().at(i));
		BOOST_CHECK_EQUAL(mut[i], individual.mut().at(i));
	}
}

/**
 * Test recombination
 */
BOOST_AUTO_TEST_CASE(recombine)
{
	int len = 10;
	ESDna<double, double> p0(len);

	vector<double> dna(len);
	for(int i = 0 ; i < len; i++)
		dna[i] = double(i);
	vector<double> mut = dna;
	ESDna<double, double> p1(dna, mut);

	//Recombine at position 0 (the offspring should be a copy of the second parent)
	auto offspring = ESDna<double, double>::recombine(p0, p1, 0);

	double error = 0;
	for(int i = 0; i < offspring.len(); i++)
		error += pow(offspring.dna().at(i) - p1.dna().at(i), 2.0);

	BOOST_CHECK_EQUAL(error, 0.0);

	//Recombine in the middle
	offspring = ESDna<double, double>::recombine(p0, p1, p0.len()/2);

	error = 0;
	for(int i = 0; i < offspring.len()/2; i++)
		error += pow(offspring.dna().at(i) - p0.dna().at(i), 2.0);
	for(int i = offspring.len()/2; i < offspring.len(); i++)
		error += pow(offspring.dna().at(i) - p1.dna().at(i), 2.0);

	BOOST_CHECK_EQUAL(error, 0.0);

}

BOOST_AUTO_TEST_SUITE_END ()
