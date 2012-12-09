#ifndef __POPULATION_H
#define __POPULATION_H

#include "Poco/Runnable.h"
#include "Poco/ThreadPool.h"
#include "Poco/RWLock.h"
#include <boost/foreach.hpp>
#include <vector>
#include <iterator>
#include <algorithm>
#include <math.h>
#include "individual.h"
#include "random.h"
#include "globals.h"
#include "worker.h"
#include <boost/smart_ptr/shared_ptr.hpp>

#define foreach         BOOST_FOREACH

typedef boost::shared_ptr<Individual> IndPtr;

using namespace std;
using namespace boost;

class Population {
  
  public:
        //vector<IndPtr> ind;
    static double *min_var;
    static double *max_var;
     Individual **ind;
    int _popsize;
    int nreal;    /* number of real paramaters */
    int nobj; /* number of objective functions */
    int ncon; /* number of constraints */
    Random *rand; /* random number generator */
    int numFronts;
    
    static int ncross;
    vector<vector<Individual *> > fronts;
    Poco::RWLock lock;
    Poco::ThreadPool pool;
    
    void simulate();
        /* inner class whose instances are intended to be executed by a thread
         * Each thread will run a number of simulations
         */
  
    Population(int popsize, int nreal, int nobj, int ncon);
    ~Population();
    void initialize();
    void print_pop(ostream& os);
    void readLimits(ifstream &ifs);
    void evaluate();
    void nondominated_sort();
    void  assignCrowdingDistance(int *indices, int start, int frontsize);
    void generateNewParent(Population *);
    void quicksortOnObjective(int *indices, int objCount, int left, int right);
    int randPartition(int *indices, int objCount, int left, int right);
    void swap( int *array, int a, int b );
    void quicksortOnCrowdingDistance(Individual **, int left, int right);
    int randPartitionDist(Individual **, int left, int right);
    void select (Population *new_pop);
    Individual* tournament(Individual *ind1, Individual *ind2);
    void sbxCrossover (Individual *parent1, Individual *parent2, Individual *child1, Individual *child2);
    void mutate_pop();
    void merge(Population *, Population *);
    void printObjectives(ostream& os);
    void qsortOnRank(int left, int right);
    int randPartitionRank(int left, int right);
};
#endif
