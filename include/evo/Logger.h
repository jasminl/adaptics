#pragma once

#include <string>
#include <tuple>

/**
 * Stores simulation results
 */
class Logger
{
public:

	/**
	 * Constructor, assigns file prefix and output folder
	 * \param folder The folder in which to store the log
	 * \param file_prefix The filename in which to store the log
	 */
	Logger(const std::string& folder, const std::string& file_name);

	/**
	 * Writes stats to file. Throws a failure exception if the file cannot be opened.
	 * \param stats The tuple of stats computed from the population.
	 */
	void dump_stats(const std::tuple<int, float, float>& stats) const;

private:
	std::string _output_file;	//!<Absolute path where file is saved
};
