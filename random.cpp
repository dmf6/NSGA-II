#include "random.h"

void Random::setSeed(int seed){
        rg.seed( seed );
    }
    
double Random::nextDouble( double lowerLimit, double upperLimit ) {
    uniform_real<> unidist(lowerLimit, upperLimit);
    variate_generator<boost::mt19937&, boost::uniform_real<> > uni(rg, unidist);
    return uni();
}

int Random::nextInt(int lowerLimit, int upperLimit) {
    uniform_int<> unidist(lowerLimit, upperLimit);
    boost::variate_generator<boost::mt19937&, boost::uniform_int<> > die(rg, unidist);
    return die();
}       
