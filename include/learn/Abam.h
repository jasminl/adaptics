#pragma once

#include <vector>

namespace learn
{

/**
 * \brief Convenience operators for learn lib
 */
class Ops
{
public:

	/**
	 * Sigmoid function (parameterized to be in the range [-1, 1])
	 * \param value Input to the sigmoid
	 * \return Sigmoid in the range [-1, 1]
	 */
	static
	double sigmoid(double value);
};

/**
 * \brief Implements an ABAM network
 */
class Abam
{

public:
	//typedefs, etc.
	enum {RANDOM, ZEROS, ONES};	//!<Mode of initialization, weight matrix

public:

	/**
	 * Allocates and initializes arrays
	 * \param dt Integration step
	 * \param lrate Learning rate
	 * \param maxd Threshold that determines when settling stops. If the absolute value difference between two settling
	 * iterations is below that threshold, settling stops.
	 * \param maxit The maximum number of settling iterations
	 * \param insize Input layer size
	 * \param outsize Output layer size
	 */
	Abam(double dt, double lrate, double maxd, int maxit,
			int insize, int outsize);

public:
	//member functions

	/**
	 * Fills in values for the weight matrix.
	 * \param Filling-in mode (RANDOM or ZEROS)
	 * \sa RANDOM, ZEROS
	 */
	void init_matrix(int mode);

	/**
	 * Computes the output layer activity by propagating forward the activity of the input layer
	 * and adding external output applied to the output layer.
	 * \param external_output External output vector applied to the output layer
	 * \return The absolute difference between the activity in the output layer before and after the iteration
	 */
	double forward(const std::vector<double>& external_output = std::vector<double>());

	/**
	 * Computes the input layer activity by propagating backward the activity of the output layer
	 * and adding external input applied to the input layer.
	 * \param external_input External input vector applied to the input layer
	 * \return The absolute difference between the activity in the input layer before and after the iteration
	 */
	double backward(const std::vector<double>& external_input = std::vector<double>());

	/**
	 * Correlation learning of weights by computing product of input and output vectors
	 * \param input Input vector
	 * \param output Output vector
	 */
	void update_weights(const std::vector<double>& input, const std::vector<double>& output);

	/**
	 * Forward and backward propagation, clamping external vectors to both the input and output layers. Useful for training.
	 * \param external_input External vector applied to the input layer
	 * \param external_output External vector applied to the output layer
	 */
	void settle(const std::vector<double>& external_input, const std::vector<double>& external_output);

	/**
	 * Computes inference by clamping an input vector and settling until either the max
	 * number of iterations has been reached or until convergence
	 * \param input Vector to clamp to the input layer
	 */
	void inference(const std::vector<double>& input);

	/**
	 * Settles to equilibrium (or up to the maximum number of iterations specified at construction)
	 * and learns the weight matrix
	 * \param input External vector applied to the input layer
	 * \param output External vector applied to the output layer
	 */
	void train(const std::vector<double>& input, const std::vector<double>&  output);

	/**
	 * Computes inference in response to an input vector and compares with an expected output vector
	 * \param external_input External vector applied to the input layer
	 * \param expected_output Expected output vector
	 * \return The absolute value error
	 */
	double test(const std::vector<double>& external_input, const std::vector<double>& expected_output);

	/**
	 * Input layer size
	 */
	int input_size() const;

	/**
	 * Output layer size
	 */
	int output_size() const;

	/**
	 * Returns a const reference to the input layer
	 */
	const std::vector<double>& input_layer();

	/**
	 * Returns a const reference to the output layer
	 */
	const std::vector<double>& output_layer();

protected:
	double _dt;						//!<Integration step
	double _lrate;					//!<Learning rate
	double _maxd;					//!<Threshold for stopping settling
	int _maxit;						//!<Max number of settling iterations

	std::vector<double> _input;		//!<Input vector
	std::vector<double> _output;	//!<Output vector
	std::vector<double> _weight;	//!<Weight matrix
};

inline
int Abam::input_size() const
{
	return (int)_input.size();
}

inline
int Abam::output_size() const
{
	return (int) _output.size();
}

inline
const std::vector<double>& Abam::input_layer()
{
	return _input;
}

inline
const std::vector<double>& Abam::output_layer()
{
	return _output;
}

}
