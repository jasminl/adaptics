#pragma once

#include <string>

/**
   Loads an XML parameter file and extracts parameters from it. Throws exceptions if unsuccessful
   \param filename The name of the file containing parameters.
*/
void load_parameters(const std::string& filename);

/**
   Reads global parameters from XML file. Throws a runtime_exception if trying to parse un unknown
   parameter.
   \param vName The name of the variable to read from the parameter file.
   \param vValue The numeric value of the parameter read.
*/
void read_parameter(const char* vName, const char* vValue);
