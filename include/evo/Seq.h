#pragma once

#include <string>
#include <list>

/**
 * \brief Low-level sequence encapsulation with bytes
 */
class Seq
{
public:

	/**
	 * Constructor
	 * \param seq_len Sequence length
	 * \param gene_len Single gene length (in number of bits)
	 */
	Seq(int seq_len, int gene_len);

	/**
	 * Destructor, frees memory allocated to data
	 */
	~Seq();

	/**
	 * Finds all occurrences of a specific tag and outputs a list of corresponding genes.
	 * \param str The string to search into
	 * \param str_len The length of the string
	 * \param gene_length The length of a gene
	 * \param tag The tag to search for which indicates the presence of a gene
	 * \return a list of genes
	 */
	static
	std::list<std::string> read_seq(const char* str, int str_len, int gene_length, const std::string& tag);

	/**
	 * Converts a sequence of int codes to a sequence of string codes
	 * \param codes The table to convert
	 * \param table_length The number of entries in the table
	 * \param code_length The length of each code
	 * \param string_table The string table where converted codes will be stored
	 */
	static
	void inplace_int2string(int** codes, int table_length, int code_length,
		char** string_table);

	/**
	 * Randomly fills the sequence with bits
	 */
	void randomize();

	/**
	 * Retrieve the (immutable) sequence data
	 * \return a const pointer to the sequence data
	 */
	const char* data() const;

private:
	int _seq_len;	//!<Sequence length
	int _gene_len;  //!<Gene length
	char* _data;	//!<Actual gene sequence
};

inline
const char* Seq::data() const
{
	return _data;
}
