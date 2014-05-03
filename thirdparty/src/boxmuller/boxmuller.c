#include "boxmuller.h"

/* The polar form of the Box-Muller transformation */
/* from Generating Gaussian Random Numbers         */
/* by Dr. Everett (Skip) F. Carter Jr.             */
/* see http://www.taygeta.com/random/gaussian.html */


double boxmuller(void)
{
  float x1, x2, w, y1, y2;

  do {
    x1 = 2.0 * drand48() - 1.0;
    x2 = 2.0 * drand48() - 1.0;
    w = x1 * x1 + x2 * x2;
  } while ( w >= 1.0 );

  w = sqrt( (-2.0 * log( w ) ) / w );
  y1 = x1 * w;
  y2 = x2 * w;

  return y1;
}

int _boxmuller_setrandomseed(void)
{
  struct timeval tv;
  unsigned int seed;

  gettimeofday(&tv, NULL);
  seed = (unsigned int) tv.tv_usec +
    (unsigned int) tv.tv_sec;

  srand48(seed);

  return 0;
}

double _boxmuller_random_double(double min, double max)
{
  return ((double) lrand48() * (max - min) / (double) RAND_MAX + min);
}

