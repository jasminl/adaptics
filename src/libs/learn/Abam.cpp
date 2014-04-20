#include "learn/Abam.h"

#include <algorithm>
#include <stdexcept>
#include <cmath>
#include <string.h>

using namespace std;

double Ops::sigmoid(double value)
{
	return -1 + 2.0 / (1 + exp(-10 * value));
}

Abam::Abam(double dt, double lrate, double maxd, int maxit,
		int insize, int outsize)
: _dt(dt), _lrate(lrate), _maxd(maxd), _maxit(maxit)
{
	_input.resize(insize, 0.0);
	_output.resize(outsize, 0.0);
	_weight.resize(insize * outsize, 0.0);
}

void Abam::init_matrix(int mode)
{
	switch (mode)
	{
	case RANDOM:
		for_each(_weight.begin(), _weight.end(), [](double& e){e = (double) rand() / RAND_MAX;});
		break;
	case ZEROS:
		_weight.assign(_weight.size(), 0);
		break;
	case ONES:
		_weight.assign(_weight.size(), 1);
		break;
	default:
		throw runtime_error("\nAbam::init_matrix: invalid mode");
	}
}

double Abam::forward(const vector<double>& external_output)
{
	auto val_before = _output;
	double sum = 0;
	for (int i = 0; i < (int) _output.size(); i++)
	{
		_output[i] -= _dt * _output[i];

		for (int j = 0; j < (int) _input.size(); j++)
			_output[i] += _dt * _weight[i * _input.size() + j] * Ops::sigmoid(_input[j]);

		if(!external_output.empty())	//Clamp external output vector
			_output[i] += external_output[i];

		sum += fabs(_output[i] - val_before[i]);
	}
	return sum;
}

double Abam::backward(const std::vector<double>&  external_input)
{
	auto val_before = _input;

	double sum = 0;

	for (int i = 0; i < _input.size(); i++)
	{
		_input[i] -= _dt * _input[i];

		for (int j = 0; j < _output.size(); j++)
			_input[i] += _dt * _weight[j * _input.size() + i] * Ops::sigmoid(_output[j]);

		if(!external_input.empty())	//Clamp external input vector
			_input[i] += external_input[i];

		sum += fabs(_input[i] - val_before[i]);
	}
	return sum;
}

void Abam::update_weights(const std::vector<double>& input, const std::vector<double>& output)
{
	for (int i = 0; i < output.size(); i++)
	{
		for (int j = 0; j < input.size(); j++)
		{
			_weight[i * input.size() + j] += _dt
					* (-_weight[i * input.size() + j]
							+ _lrate * Ops::sigmoid(output[i])
									* Ops::sigmoid(input[j]));
		}
	}
}

void Abam::settle(const vector<double>& external_input, const vector<double>& external_output)
{
	int iter = _maxit;
	while (iter--)
	{	//Note the backward pass is evaluated only if the forward is below convergence threshold
		if(forward(external_output) < _maxd && backward(external_input) < _maxd)
			break;
	}
}

void Abam::inference(const vector<double>& input)
{
	int iter = _maxit;
	while (iter--)
	{
		if (forward() < _maxd && backward(input) < _maxd)
			break;
	}
}

void Abam::train(const vector<double>&  input, const vector<double>&  output)
{
	_output.assign(_output.size(), 0);
	_input = input;

	settle(input, output);
	update_weights(_input, _output);
}

double Abam::test(const std::vector<double>& external_input, const std::vector<double>& expected_output)
{
	_output.assign(_output.size(), 0);
	_input = external_input;

	inference(external_input);

	double error = 0;
	for(int i = 0; i < (int)_output.size(); i++)
		error += fabs(Ops::sigmoid(_output[i]) - Ops::sigmoid(expected_output[i]));

	return error;
}
