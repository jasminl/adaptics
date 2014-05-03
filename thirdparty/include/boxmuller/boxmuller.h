#ifndef _BOXMULLER_H
#define _BOXMULLER_H

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <sys/time.h>

double boxmuller(void);
int _boxmuller_setrandomseed(void);
double _boxmuller_random_double(double min, double max);


#endif
