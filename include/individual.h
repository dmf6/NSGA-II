#ifndef __INDIVIDUAL_H
#define __INDIVIDUAL_H

#include <boost/smart_ptr.hpp>
#include <vector>
#include "random.h"

using namespace std;
using namespace boost;

class Individual {
  public:
    double *objs; /* fitness functions */
    double *vars; // stores params
    double *xreal; //stores zfstats
    double *constr;
    double crowding_distance;
    int rank;
    int numvars;
    int nobj;
    int ncon;
    int ndom;    // domination count
    vector<Individual *> ss; // set of solutions dominated by this solution
    
    Random *rand;
    double constr_violation;
    int index;
    
    Individual(int nreal, int nobj, int ncon);
    ~Individual();
    void initialize(Random *rand, double *, double *);
    void evaluate();
    int checkDominance(Individual *);
    void mutate(double *, double *);
    void simulate(int id);
    void copy(Individual *);
};
#endif
