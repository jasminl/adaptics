#pragma once

extern "C"
{
	#include "boxmuller.h"
}

#include <vector>

namespace evo
{
	/**
	 * Base class for ES sequences
	 */
	template <typename T, typename G>
	class DnaCont
	{
	public:
		//typedefs, etc.
		static const int DEF_REC_PT = -1;		//!<Indicate to automatically (randomly) select a recombination point

	public:

		/**
		 * Default constructor
		 */
		DnaCont(){}

		/**
		 * Allocates dna and mutation operator sequences
		 * \param seq_length Length of DNA sequence to allocate
		 */
		DnaCont(int seq_length)
		{
			_dna.resize(seq_length, T());
			_mut.resize(seq_length, G());
		}

		/**
		 * Allocates dna and mutation operator and sets their values
		 * \param dna DNA sequence to copy
		 * \param mut Mutation variable sequence to copy
		 */
		DnaCont(std::vector<T> dna, std::vector<G> mut)
		{
			//Copy-and-swap idiom
			using std::swap;
			std::swap(_dna, dna);
			std::swap(_mut, mut);
		}

		/**
		 * \return Length of the DNA sequence
		 */
		int len() const
		{
			return (int)_dna.size();
		}

		/**
		 * \return A const reference to the dna
		 */
		const std::vector<T>& dna() const
		{
			return _dna;
		}

		/**
		 * \return A const reference to the mutation vector
		 */
		const std::vector<G>& mut() const
		{
			return _mut;
		}

	protected:
		std::vector<T> _dna;	//!<The genotype
		std::vector<G> _mut;	//!<Characteristics of the mutation operator
	};

	/**
	 * Evolution strategies template specialization for double-precision sequences
	 * with double-precision Gaussian spread parameters.
	 */
	template <typename T, typename G>
	class ESDna
	: public DnaCont<T, G>
	{
	public:

		/**
		 * Allocates sequence
		 * \param seq_length Length of DNA sequence to allocate
		 * \sa DnaCont
		 */
		ESDna(int seq_length)
		: DnaCont<T, G>(seq_length)
		  {}

		ESDna(std::vector<T> dna, std::vector<G> mut)
		: DnaCont<T, G>(dna, mut)
		{}


		template <typename K, typename H>
		void mutate_in_place(std::vector<K>& dest, std::vector<H>& src)
		{
			//Multiples (0, 1)-distributed Gaussian number by item-specific spread
			for(int i = 0; i < (int) dest.size(); i++)
				dest[i] += src[i] * boxmuller();
		}

		/**
		 * Apply Evolution Strategies mutation operator where the mutation is given by a
		 * Gaussian-distributed random number
		 */
		void mutate()
		{
			mutate_in_place(this->_dna, this->_mut);
		}

		/**
		 * Mutates the mutation operator
		 */
		void mutate_operator()
		{
			mutate_in_place(this->_mut, this->_mut);
		}


		template<typename K>
		static
		std::vector<K> recombine(const std::vector<K>& first, const std::vector<K>& second, int rec_pt)
		{
			if(first.size() != second.size())
				throw std::runtime_error("ESDna<double, double>::recombine_in_place: must be of equal lengths");

			int point = rec_pt;
			if(rec_pt == DnaCont<T,G>::DEF_REC_PT)
				point = rand() % (int)first.size();

			auto out = first;
			memcpy(&out[point], &second[point], sizeof(T) * (second.size() - point));
			return out;
		}

		static
		ESDna<T, G> recombine(const ESDna<T, G>& parent0, const ESDna<T, G>& parent1,
					int rec_pt = DnaCont<T,G>::DEF_REC_PT)
		{
			return {recombine(parent0.dna(), parent1.dna(), rec_pt), std::vector<G>(parent1.len())};
		}

		void recombine_operator(const ESDna<T, G>& parent0, const ESDna<T, G>& parent1,
					int rec_pt = DnaCont<T,G>::DEF_REC_PT)
		{
			if(parent0.len() != parent1.len())
				throw std::runtime_error("ESDna<double, double>::recombine: parents must be of equal lengths");

		}
	};
}
