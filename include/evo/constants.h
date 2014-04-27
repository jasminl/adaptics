#pragma once

namespace evo
{

/*Maximum number of iterations before construction of different individuals. If individuals and blocks are not different after 
MAX_ITERATIONS trials, the program ends. */
#define MAX_ITERATIONS 50	

// Means all individuals are different
#define DIFFERENT -1

// Means the construction does not converge
#define NOTCONVERGE -2

//Maximum fitness value that can be achieved by an individual
#define OPTIMUM	31

/////////////////////////////////////////////////////////////////////////

//maximum iteration for different blocks formation
#define MAX_ITERATIONS2 50

//means blocks are all different
#define DIFFERENT2 -1

//means it is not possible to create all different blocks
#define NOTCONVERGE2 -2

//means no building blocks were found
#define NOTFOUND -3

//to continue simulating
#define CONTINUE 10210209

#define OPTIMAL_INDIV 987987986

}
