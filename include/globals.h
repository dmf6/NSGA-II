#ifndef __GLOBALS_H
#define __GLOBALS_H

#include <cstdlib>
#include <fstream>
#include <sstream>
#include <iostream>

/*popsize must be a multiple of 4 */
#define POP_SIZE 400
#define NUM_PARS 1
#define NUM_OBJS 2
#define NUM_GEN 100
#define NUM_CON 0

#define INF 1.0e14                              
#define EPS 1.0e-14
/* crossover distribution index between 5 and 20 */
#define ETA_C 20
/* mutation distribution index between 5 and 50 */
#define ETA_M 20
#define P_CROSS 0.9
#define P_MUT 0.5

#define NUM_THREADS 16
/* 16 threads working on a poopulation of size 400 means that each
 * thread will run 25 simulations */
#define NUM_SIMS 50
#define PI 3.141592653589793

void costfunc(double *, double *, double *);

#endif

