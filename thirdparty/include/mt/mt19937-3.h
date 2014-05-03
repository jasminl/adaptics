/*
 * Header file generated to encapsulate Makoto Matsumoto and Takuji Nishimura's
 * mt19937-3.c source file so as to compile it into a third party library.
 *
 * @author: Jasmin Leveille
 * @email: jalev51@gmail.com
 */

/*
 * Initializing the array with a seed
 * \param seed The seed to initialize the random generator to.
 */
void
sgenrand(unsigned long seed);

void
lsgenrand(unsigned long seed_array[]);

/**
 * Generating reals
 * \return A double-precision random number
 */
double
genrand();
