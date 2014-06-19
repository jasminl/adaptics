#pragma once

#include <map>
#include <string>

namespace evo
{

static const int gray_undefined = -1; //!<For a string that has no corresponding code

/**
 * Returns the next position where to insert a bit. Throws a runtime_error exception in case of an invalid
 * return value, which could result from bits not having been initialized with UNALLOCATED.
 * \param line Sequence to modify
 * \param code_length Length of the sequence
 * \return The position where to insert the bit
 */
int find_cur_pos(int* line, int code_length);

/**
 * Copies a subset of bits from one sequence to another.
 * \param dest Destination sequence where bits will be copied to.
 * \param src Source sequence where bits will be copied from.
 * \param start Starting position in the sequence where to start copying bits
 * \param code_length Number of bits in a sequence
 */
void inplace_relocate(int* dest, int* src, int start, int code_length);

/**
 * Creates a Gray code table.
 * \param codes The resulting array of Gray code sequences
 * \param table_length The number of sequences in the table
 * \param code_length The length of each sequence
 * \return A map of gray code strings and corresponding integers
 */
std::map<std::string, int> create_gray(int** codes, int table_length, int code_length);

/**
 * Prints the Gray code to file
 * \param file_path Path of file where to print the Gray code
 * \param codes Gray code table to print out
 * \param table_length The number of sequences in the table
 * \param code_length The length of each Gray code sequence
 * \sa create_gray
 */
void print_table(const std::string& file_path, int** codes, int table_length, int code_length);

/**
 * Match a graycode
 * \param target The code to locate
 * \param table The table of codes
 */
int match_gray(const std::string& target, const std::map<std::string, int>& table);

}
