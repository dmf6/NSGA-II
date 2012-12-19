#include "population.h"
#include <map>

int Population::ncross = 0;

 bool compareDistance(const Individual *a, const Individual *b) {
    return (a->crowding_distance > b->crowding_distance);
}

/* inline constructor */
Population::Population(int popsize, int nreal, int nobj, int ncon) : _popsize(popsize) {
        /* need to seed from the random number generator for the entire population */
    rand = new Random();
    rand->setSeed(static_cast<unsigned int>(time(0)));
    this->nreal = nreal;
    this->nobj = nobj;
     ind = new Individual*[_popsize];
        /* since we use boost shared pointers no need to worry about
         * freeing up memory as it is automatic */
    for (int i = 0; i < _popsize; i++) {
        // IndPtr indPtr(new Individual(nreal, nobj, ncon));
        // ind.push_back(indPtr);
        
        ind[i] = new Individual(nreal, nobj, ncon);
        //     /* let each individual have an identification number We
        //      *  need this for sorting the indices for each front based
        //      *  on each objective
        //      */
         ind[i]->index = i;
    }
    ncross = 0; //static member so this value will apply for all class instances
    this->ncon = ncon;
    numFronts = 0;
}

Population::~Population() {
    for (int j = 0; j < _popsize; j++) {
        delete ind[j];
    }
    delete [] ind;
/* no need to delete static arrays since they are not allocated on the
 * stack or heap */
    delete rand;
}

/* initialize all of the individuals already created in constructor
 * with random parameters */ 
void Population::initialize() {
    for (int i = 0; i < _popsize; i++) {
        ind[i]->initialize(rand, min_var, max_var);     
    }
}

void Population::print_pop(ostream& os) {
    os  << "Population size = " << _popsize << "\t" << " Number of Parameters being optimized = " << nreal << endl;
    
    for (int i = 0; i < _popsize; i++) {
        for(int j = 0; j < nreal; j++) {
            os << " " << ind[i]->vars[j] << " ";
        }
        os << "\n";   
    }
    os << endl;
    
}

void Population::printObjectives(ostream& os) {
    for (int i = 0; i < _popsize; i++) {
        for(int j = 0; j < nobj; j++) {
            os << " " << ind[i]->objs[j] << " ";
        }
        os << "\n";   
    }
    os << endl;
}

void Population::evaluate() {
    for (int j = 0; j < _popsize; j++) {
            // each individual in the popultion will be evaluated against each objective
        Individual *current_ind = ind[j];
        current_ind->evaluate();
    }
}

// /*need to store a array of indicies for the population this is the
//  *  array that we sort using randomized qsort and then we use the
//  *  elements of this array to index the individuals in the current
//  *  front to calculate crowding distance */

void Population::nondominated_sort() {
        //clear fronts list for the current generation
    fronts.clear();
    vector< Individual *> front;    /* first front */
    int *indices = new int[_popsize];
    int icount = 0;
    
     int flag = 0;
    
    for (int p = 0; p < _popsize; p++) {
        Individual *ind_p = ind[p];
        for (int q = 0; q < _popsize; q++) {
            Individual *ind_q = ind[q];
                // check to see if p dominates q
            flag = ind_p->checkDominance(ind_q);      
            if (flag==1) {
                ind_p->ss.push_back(ind_q);
            } else if (flag==-1){
                    //update domination count for individual 
                ind_p->ndom++;
            }    
        }     
        if (ind_p->ndom == 0) {
                /* If there are no individuals that dominate ind_p then it belongs to the
                 *  first front. Each individual has a member variable called rank
                 */
                //cout << "No individuals dominate individual " << ind_p->index << endl;
            
            ind_p->rank = 1;
                //store individual in first front because it is 'non-dominated' by all solutions
            front.push_back(ind_p);
                /* store index of element too, first front element indices come first obviously */
            indices[icount++] = ind_p->index;
        }
    }
    /* store first front */
    fronts.push_back(front);
    
    assignCrowdingDistance(indices, 0, front.size());
    
    int start = front.size();
    
    numFronts = 1; //representing the first front
    vector<Individual *> q; /* used to store members of the next front */
    while(!front.empty()) {
            //cout << "Creating a new front..." << endl;
        q.clear();
            /* iterate through all members of the first front */
        foreach(Individual *ind, front) {
            vector<Individual *> s = ind->ss;
                /* iteratre through the set of solutions dominated by it */
            foreach(Individual *qi, s) {
                qi->ndom -=1;
               
                if(qi->ndom == 0) {
                        //cout << "domination count for individual with idx " << qi->index  << " is " << qi->ndom << endl;
                    qi->rank = numFronts + 1;
                    q.push_back(qi);
                        //cout << "The number of individuals in q is " << q.size() << "\n";
                        /* store index of individual whose domination
                           count has reached 0. This individual will
                           not be visited again */
                    indices[icount++] = qi->index;
                }
            }
            front.pop_back();
        }
      
        if (!q.empty()) {
            fronts.push_back(q);
                /* assign crowding distance in current front */
            assignCrowdingDistance(indices, start, q.size());
            start = start + q.size();
            numFronts++; //increment i to identify the next front
        }
        front.swap(q);
    }
    qsortOnRank(0, _popsize-1);
     
    cout << "FINISHED NON-DOMINATED SORTING WITH " << numFronts << "sets\n";
    delete [] indices;
}

/*randomized quicksort to sort front based on objection objcount */
void Population::qsortOnRank(int left, int right) {
        //initialize pivot index
    int pivot;
    if (left < right) {
        pivot = randPartitionRank(left, right);
        qsortOnRank(left, pivot-1);
        qsortOnRank(pivot+1, right);
    }
}

/* we don't pass the population ind because it is a member of this class */
int Population::randPartitionRank(int left, int right) {
    int pivot;
    int randIdx = rand->nextInt(left, right);

    std::swap(ind[right], ind[randIdx]);

    pivot = ind[right]->rank;

     int k = left - 1;
     for (int j = left; j < right; j++) {
         if(ind[j]->rank <=pivot) {
              k++;
              std::swap(ind[k], ind[j]);
         }
     }
     std::swap(ind[k+1], ind[right]);
    return k+1;
}


/* notice that the population members aren't actually being switched
 * just the indices array and this is used to index the population
 * when calculating the crowding distance */
void Population::assignCrowdingDistance(int *indices, int start, int frontsize){
    if (frontsize == 1){
        ind[indices[0]]->crowding_distance = INF; 
        return;
    }
    if (frontsize==2) {
        ind[indices[0]]->crowding_distance = INF;
        ind[indices[1]]->crowding_distance = INF;
        return;
    }
        /* for each objective call randqsort to sort the set in worse
         * order of current obejctive (return sorted indices array) */
    
    for (int i = 0; i < nobj; i++) {
            //using randomized quicksort, sort the front indices array
            // in order of objective i use the sorted indices array to
            // index the population when calculating crowding distance
            // for all members in front
  
        quicksortOnObjective(indices, i, 0, frontsize-1);

            /* assign large distance to boundary solutions */
        ind[indices[0]]->crowding_distance = INF;
        ind[indices[frontsize-1]]->crowding_distance = INF;
        double fmax = ind[indices[frontsize-1]]->objs[i];
        double fmin = ind[indices[0]]->objs[i];
        
                /* all other solutions 2-(n-1) */
        double nextObj, previousObj;
        for (int j = 1; j < frontsize-1; j++) {
            nextObj = ind[indices[j+1]]->objs[i];
            previousObj =  ind[indices[j-1]]->objs[i];
            if ((fmax - fmin) == 0) {
                ind[indices[j]]->crowding_distance = INF;
            }
            else {
                ind[indices[j]]->crowding_distance += (nextObj - previousObj)/( fmax - fmin);
            }    
        }
    }
    for (int k = 0; k < frontsize; k++) {
        if (ind[indices[k]]->crowding_distance != INF) {
            ind[indices[k]]->crowding_distance = (ind[indices[k]]->crowding_distance)/nobj;
        }
    }
    
}

/*randomized quicksort to sort front based on objection objcount */
void Population::quicksortOnObjective(int *indices, int objCount, int left, int right) {
        //initialize pivot index
    int pivot;
    if (left < right) {
            //pass objective number
        pivot = randPartition(indices, objCount, left, right);
        quicksortOnObjective(indices, objCount, left, pivot-1);
        quicksortOnObjective(indices, objCount, pivot+1, right);
    }
}

/* we don't pass the population ind because it is a member of this class */
int Population::randPartition(int *indices, int objCount, int left, int right) {
    double pivot;
    int randIdx = rand->nextInt(left, right);
    swap(indices, right, randIdx);
        /* use random index as pivot and use indices to index the
         * population and objCount to index objective array*/
    pivot = ind[indices[right]]->objs[objCount];
    int k = left - 1;
    for (int j = left; j < right; j++) {
        if(ind[indices[j]]->objs[objCount] <=pivot) {
            k++;
            swap(indices, k, j);
        }
    }
    swap(indices, k+1, right);
    return k+1;
}

void Population::swap( int *array, int a, int b ) {
    int temp = array[ a ];
    array[ a ] = array[ b ];
    array[ b ] = temp;
}


void Population::select (Population *new_pop) {
        /* need to initialize before using */
    Individual *parent11 = NULL;
    Individual *parent12 = NULL;
    
    vector<int> a1(_popsize);
    // vector<int> a2(_popsize);
    
    for (int i=0; i<_popsize; i++) {
        a1[i] = i;
    }
        /* shuffle to ensure the same individual doesn't mate with
         * itself. There is still a chance but the higher the
         * population size, the more unlikely this is. */
    random_shuffle ( a1.begin(), a1.end() );
    // random_shuffle ( a2.begin(), a2.end() );
    double randIdx1, randIdx2 = 0;
    
        /* let the parents fight it out to get it on but don't let them mate with themselves */
    int count = 0;

        // keep a temp vector of individuals chosen to be parents
    map<Individual *, Individual *> parents;

     randIdx1 = rand->nextInt(0, a1.size()-1);
     randIdx2 = rand->nextInt(0, a1.size()-1);
        
     parent11 = tournament(ind[a1[randIdx1]], ind[a1[randIdx2]]);
     parent12 = tournament(ind[a1[randIdx1]], ind[a1[randIdx2]]);
     while (parent11 == parent12) {
         cout << "HELLO2\n";
         randIdx1 = rand->nextInt(0, a1.size()-1);
         randIdx2 = rand->nextInt(0, a1.size()-1);
         parent12 = tournament(ind[a1[randIdx1]], ind[a1[randIdx2]]);
     }
         /* remember these two parents */
     parents[parent11] = parent12;
         /*crossover ONCE */
     sbxCrossover (parent11, parent12, new_pop->ind[count], new_pop->ind[count+1]);
     count = count + 2;

     cout << "PARENTS ADDED " << parent11->index << " " << parent12->index << "\n";
     cout << count << ":RESPECTIVE RANK AND CROWDING DISTANCE " << parent11->rank << " " << parent11->crowding_distance << " " <<  parent12->rank << " " << parent12->crowding_distance << "\n";
     while (count < _popsize) {
         cout << "HELLO1\n";
         
         randIdx1 = rand->nextInt(0, a1.size()-1);
         randIdx2 = rand->nextInt(0, a1.size()-1);
        
         parent11 = tournament(ind[a1[randIdx1]], ind[a1[randIdx2]]);
         parent12 = tournament(ind[a1[randIdx1]], ind[a1[randIdx2]]);
         cout << parent11->index << " " << parent12->index << "\n";
             /* we want distinct parents */
       
             /* if parent 1 is not already in the list then add */
         if (parents.find(parent11) == parents.end()) {
             while ((parents.find(parent12) != parents.end() || (parent11 == parent12))) {
                 randIdx1 = rand->nextInt(0, a1.size()-1);
                 randIdx2 = rand->nextInt(0, a1.size()-1);
                 parent12 = tournament(ind[a1[randIdx1]], ind[a1[randIdx2]]);
             }
             cout << count <<":PARENTS ADDED " << parent11->index << " " << parent12->index << "\n";
             cout << "RESPECTIVE RANK AND CROWDING DISTANCE " << parent11->rank << " " << parent11->crowding_distance << " " <<  parent12->rank << " " << parent12->crowding_distance << "\n";
             
             cout << "HELLO3\n";
             parents[parent11] = parent12;
             sbxCrossover (parent11, parent12, new_pop->ind[count], new_pop->ind[count+1]);
             count = count + 2;
         }
         else {
             while ((parents.find(parent11)->first) == parent11) {
                 cout << "HELLO4\n";
                 cout << parent11->index << " " << parent12->index << "\n";
                 randIdx1 = rand->nextInt(0, a1.size()-1);
                 randIdx2 = rand->nextInt(0, a1.size()-1);
                 parent11 = tournament(ind[a1[randIdx1]], ind[a1[randIdx2]]);
             }
             
            
            /* if parent11 has already fought with parent12 then choose another parent12*/
             while (((parents.find(parent11)->second) == parent12) || ((parents.find(parent11)->first) == parent12) || (parent11 == parent12)) {
                    // keep parent11 but change parent12 because we
                    // have see parent11 but we just want to let it
                    // fight with another individual, notice that it
                    // could still dominate another individual and so
                    // will get passed on twice */
                randIdx2 = rand->nextInt(0, a1.size()-1);
                parent12 = tournament(parent11, ind[a1[randIdx2]]);
                cout << parent11->index << " " << parent12->index << "\n";
            }
             
            parents[parent11] = parent12;
                /* If we have reached this point then these two
                 * parents have never mated before and they are not
                 * the same */
            sbxCrossover (parent11, parent12, new_pop->ind[count], new_pop->ind[count+1]);
            count = count + 2;
            cout << "PARENTS ADDED! " << parent11->index << " " << parent12->index << "\n";
            cout << count << ":RESPECTIVE RANK AND CROWDING DISTANCE " << parent11->rank << " " << parent11->crowding_distance << " " <<  parent12->rank << " " << parent12->crowding_distance << "\n";
         }
     }
}


Individual* Population::tournament (Individual *ind1, Individual *ind2) {
    double crowd1 = ind1->crowding_distance;
    double crowd2 = ind2->crowding_distance;
    int rank1 = ind1->rank;
    int rank2 = ind2->rank;

    if (rank1 < rank2){
        return (ind1);
    }
    else if (rank2 < rank1) {
        return (ind2);
    }
    else if (crowd1 > crowd2) {
        return(ind1);
    }
    else if (crowd2 > crowd1) {
        return(ind2);
    }
    else if (rand->nextDouble(0, 1) <= 0.5) {
        return(ind1);
    }
    else {
        return(ind2);
    }
}

/* Routine for real variable SBX crossover */
void Population::sbxCrossover (Individual *parent1, Individual *parent2, Individual *child1, Individual *child2) {
    double random;
    double y1, y2, yl, yu;
    double c1, c2;
    double alpha, beta, betaq;
    
    if (rand->nextDouble(0, 1)<=P_CROSS) {
        for (int i=0; i<NUM_PARS; i++) {
            if (rand->nextDouble(0, 1) <= 0.5) {
                 if (fabs(parent1->vars[i]-parent2->vars[i]) > EPS) {
                     if (parent1->vars[i] < parent2->vars[i]) {
                         y1 = parent1->vars[i];
                         y2 = parent2->vars[i];
                     }
                     else {
                         y1 = parent2->vars[i];
                         y2 = parent1->vars[i];
                     }
                     yl = min_var[i];
                     yu = max_var[i];

                     random = rand->nextDouble(0, 1);
                     beta = 1.0 + (2.0*(y1-yl)/(y2-y1));
                     alpha = 2.0 - pow(beta,-(ETA_C+1.0));
                     if (random<= (1.0/alpha)) {
                         betaq = pow ((random*alpha),(1.0/(ETA_C+1.0)));
                     }
                     else {
                         betaq = pow ((1.0/(2.0 - random*alpha)),(1.0/(ETA_C+1.0))); 
                     }
                     c1 = 0.5*((y1+y2)-betaq*(y2-y1));                
                     beta = 1.0 + (2.0*(yu-y2)/(y2-y1));
                     alpha = 2.0 - pow(beta,-(ETA_C+1.0));
                     if (random <= (1.0/alpha)) {
                         betaq = pow ((random*alpha),(1.0/(ETA_C+1.0)));
                     }
                     else {
                         betaq = pow ((1.0/(2.0 - random*alpha)),(1.0/(ETA_C+1.0)));
                     }
                     c2 = 0.5*((y1+y2)+betaq*(y2-y1));
                         /* Make sure that the generated element is
                          * within the specified decision space else
                          * set it to the appropriate extrema */
                     if (c1<yl) {
                         c1=yl;
                    }
                     else if (c1>yu) {
                         c1=yu;
                     }
                     if (c2<yl){
                         c2=yl;
                     }
                     else if (c2>yu) {
                        c2=yu;
                    }
                     
                     if (rand->nextDouble(0, 1) <= 0.5) {
                        child1->vars[i] = c2;
                        child2->vars[i] = c1;
                    }
                    else {
                        child1->vars[i] = c1;
                        child2->vars[i] = c2;
                    }
                 }
                 
                 else {
                     child1->vars[i] = parent1->vars[i];
                     child2->vars[i] = parent2->vars[i];
                 }
            }
            else {
                child1->vars[i] = parent1->vars[i];
                child2->vars[i] = parent2->vars[i];
            }
        }
    }
    else {
        for (int j=0; j < NUM_PARS; j++) {
                child1->vars[j] = parent1->vars[j];
                child2->vars[j] = parent2->vars[j];
        }
    }
}

// void Population::simulate() {
//         /* Parallel implementation */
//       vector<Poco::Runnable *> runnables;
//         for (int i = 0, j = 0; i < NUM_THREADS; i++, j+=NUM_SIMS) {
//             runnables.push_back(new Worker(i, &lock, j, j+NUM_SIMS, this));
//         }
//         foreach( Poco::Runnable *run, runnables) {
//             pool.start(*run);
//         }
//         pool.joinAll();
// }

void Population::merge(Population *pop1, Population *pop2) {
        /* copy parent population into positions 1-popsize */
        /* remember mixed pop is twice the size as parent and child pop */
    for (int i = 0; i < pop1->_popsize; i++) {
    //     /* copy individuals into new mixed population */
        ind[i]->copy(pop1->ind[i]);
    }     
        /*copy child population into positions popsize-2*popsize */
    for (int i = 0, j = pop2->_popsize; j < _popsize; i++, j++) {
        ind[j]->copy( pop2->ind[i]);
    }
}

void Population::mutate_pop() {
    for (int i=0; i<_popsize; i++) {
            /* need to pass in lower and upper limit to polynomial mutation function */
        ind[i]->mutate(rand, min_var, max_var);
    }
}
    
void Population::generateNewParent(Population *parent) {
    int frontCounter = 1;
    int j = 0;

    int currentSize = 0;

    while(frontCounter < numFronts) {
        vector<Individual *> front;
        front.clear();
        while(j < parent->_popsize) {
            if(ind[j]->rank == frontCounter) {
                front.push_back(ind[j]);
                j++;
            }
            else {
                frontCounter++;
                break;
            }
        }
        if (currentSize + (int) front.size() < parent->_popsize) {
                //cout << "Front size is " << currentSize << "\n";
            for (int p = currentSize, q =0 ; q < (int) front.size(); q++, p++) {
                parent->ind[p]->copy(front[q]);
            }
            
            currentSize += front.size(); 
        }
        else {
            sort(front.begin(), front.end(), compareDistance);
    
            for (int m =0, n = currentSize; n < parent->_popsize; m++, n++) {
                parent->ind[n]->copy(front[m]);
            }
            break;
        }
        
    }
    
}


 void Population::readLimits(ifstream & ifs) {
    string line;
    double *buff = new double[2];
     if (ifs.is_open()) {  
         int i = 0;
         while (ifs.good()) {  
             getline( ifs, line);
             stringstream sstr( line );
             double value;
             int j = 0;
             
             while( sstr >> value) {
                 buff[j] =  value;
                 j++;          
             }
                 //cout << buff[0] << "\t" << buff[1] << endl;
            
             min_var[i] = buff[0];
             max_var[i] = buff[1];
             i++;
         }
         delete [] buff;
     }
     
     else {
         cout << "Unable to open limits file" << endl;
     }
}

