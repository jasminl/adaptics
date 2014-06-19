/**
 * \file Utils.h
 * \brief Misc functions to print matrices
 */
#pragma once

#include <string>
#include <fstream>
#include <iomanip>
#include <exception>

namespace visio
{
/** @name Miscellaneous functions
*/
//@{

/**
 * Dump matrix to file
 * @param fname File name to save into
 * @param matrix Matrix to save
 * @param vsize Number of rows in matrix
 * @param hsize Number of columns in matrix
 */
template <typename T>
void print_matrix(const std::string& fname, T** matrix, unsigned int vsize, unsigned int hsize)
{
	std::ofstream ofile(fname);
	if(ofile.bad())
		throw std::ios_base::failure("print_matrix: unable to open file for writing " + fname);

	for(unsigned int i = 0; i < vsize; i++)
	{
		for(unsigned int j = 0; j < hsize; j++)
			ofile<<matrix[i][j]<<" ";
		ofile<<std::endl;
	}
}

template <typename T>
void printMatrix(const std::string& fname,T** matrix, unsigned int trueVsize,
		unsigned int trueHsize, unsigned int vsize)
{
	std::ofstream ofile(fname);

	for(unsigned int i = 0; i < vsize; i++)
	{
        for(unsigned int j = 0; j < trueVsize; j++)
		{
			for(unsigned int k = 0; k < trueHsize; k++)
				ofile<<std::setprecision(3)<<matrix[i][j * trueHsize + k]<<" ";
			ofile<<std::endl;
		}
		ofile<<std::endl;
	}
}
//@}
}
