#include <iostream>
#include <fstream>
#include <cstdlib> 
#include "individual.h"
#include "globals.h"

using namespace std;

Individual::Individual(int nreal, int nobj, int ncon) {
    this->nobj = nobj;
    numvars = nreal;
    vars = new double[nreal];
    objs = new double[nobj];
    xreal = new double[nobj];
    constr = new double[ncon];
    memset(vars, 0.0, nreal);
    memset(objs, 0.0, nobj);
    memset(xreal, 0.0, nobj);
    memset(constr, 0.0, ncon);
   
    this->ncon = ncon;
    crowding_distance = 0;
    constr_violation = 0;
    rank = 0;
    ndom = 0;
}

Individual::~Individual() {
    delete [] vars;
    delete [] objs;
    delete [] xreal;
    delete [] constr;
}

void Individual::initialize(Random *rand, double *min, double *max) {
    double ran;
    for (int j=0; j<numvars; j++) {
        ran= rand->nextDouble(min[j], max[j]);
        vars[j] = ran;
    }
}

void Individual::evaluate() {
        // vars is the parameter array
        //xreal is the objective measurement
        //costfunc(&xreal[0], &objs[0], &constr[0]);
    costfunc(vars, objs, constr);
        //cout << "Objective values are " << objs[0] << " " << objs[1] << endl;
    constr_violation = 0.0;
    for (int j=0; j<ncon; j++) {
        if (constr[j]<0.0) {
            cout << "There is a constraint violation " << endl;
            constr_violation += constr[j];
        }
    }
}

/* Comparator function to check for dominance
 *
 * this instance dominates 'another' instance if
 * (1) the solution is no worse than 'another' solution in
 * all objectives (1, ..., NOBJ)
 * OR
 * This solution is better than another solution in at least one objective 
 */
int Individual::checkDominance(Individual *another) {
    int flag1, flag2;
    flag1 = 0; flag2 = 0;

    if (constr_violation<0 && another->constr_violation<0) {
        if (constr_violation > another->constr_violation) {
                    /* this and another are both infeasible but this has a smaller overall constraint violation */
            return (1);
        }
        else {
            if (constr_violation < another->constr_violation) {
                 /* this and another are both infeasible but another has a smaller overall constraint violation */
                return (-1);
            }
            else {
                return (0);
            }
        }
    }
    else {
        if (constr_violation < 0 && another->constr_violation == 0) {
            return (-1);
        }
        else {
            if (constr_violation == 0 &&  another->constr_violation < 0) {
                return (1);
            }
            else {
                    /* if we get to this point then there are no
                     * constraint violations and we both solutions'
                     * objective vectors are compared */
                for (int i=0; i < NUM_OBJS; i++) {
                    if (objs[i] < another->objs[i]) {
                        flag1 = 1;
                            // cout << "This Individual dominates another individual" << endl;
                    }
                    else if (objs[i] > another->objs[i]){
                        flag2 = 1;
                        // cout << "This Individual is dominated by another individual" << endl;
                    }       
                }
            }
            
            if (flag1==1 && flag2==0) {
                return (1);
            }
            else {
                if (flag1==0 && flag2==1) {
                    return (-1);
                }
                else {
                    return (0);
                }
            }
        }
    }
}


// /* Routine for real polynomial mutation of an individual */
void Individual::mutate(Random *rand, double *min_var, double *max_var) {
    double rnd, delta1, delta2, mut_pow, deltaq;
    double y, yl, yu, val, xy;
    
    for (int j=0; j<NUM_PARS; j++) {
        rnd= rand->nextDouble(0, 1);
        if ( rnd <= P_MUT) {
            y = vars[j];
            yl = min_var[j];
            yu = max_var[j];
            delta1 = (y-yl)/(yu-yl);
            delta2 = (yu-y)/(yu-yl);
            rnd = rand->nextDouble(0, 1);
            mut_pow = 1.0/(ETA_M+1.0);
            if (rnd <= 0.5) {
                xy = 1.0-delta1;
                val = 2.0*rnd+(1.0-2.0*rnd)*(pow(xy,(ETA_M+1.0)));
                deltaq =  pow(val,mut_pow) - 1.0;
            }
            else {
                xy = 1.0-delta2;
                val = 2.0*(1.0-rnd)+2.0*(rnd-0.5)*(pow(xy,(ETA_M+1.0)));
                deltaq = 1.0 - (pow(val,mut_pow));
            }
            y = y + deltaq*(yu-yl);
            if (y<yl) {
                y = yl;
            }
            if (y>yu) {
                y = yu;
            }
            vars[j] = y;
        }
    }
}

// /* id is the thread number */
void Individual::simulate(int id) {
        //convert double array into parameter string
        //could use either stringstream or boost::lexical_cast (faster)

    string prog = "./neuralsim -i";
    string sid = boost::lexical_cast<string>(id);
    string params = " -p ";
    string par;
            
    for (int j = 0; j < 21; j++) {
        par = boost::lexical_cast<string>(vars[j]);
        params += " " + par;
    }
            
    prog += " " + sid;
    prog += params;
        
    cout << prog << endl;
    
    boost::scoped_array<char> writable(new char[prog.size() + 1]);
    std::copy(prog.begin(), prog.end(), writable.get());
    writable[prog.size()] = '\0'; // don't forget the terminating 0
    int i;
            
    i = system(writable.get());
    if (i == 0) {
        cout << "Call to neuralsim successful" << endl;
    }
        /* call zfgenerator */
    string prog2 = "./zfgenerator -f ";
    string name = "SC";
    prog2 += name;
    prog2 += sid;
    prog2 += ".asc";

    cout << prog2 << endl;
    
    boost::scoped_array<char> writable2(new char[prog2.size() + 1]);
    std::copy(prog2.begin(), prog2.end(), writable2.get());
    writable2[prog2.size()] = '\0'; // don't forget the terminating 0
    i = system(writable2.get());
    if (i == 0) {
        cout << "Call to zfgenerator successful" << endl;
    }

    string myfile = "SC";
    myfile += sid;
    myfile += "_zfstats.dat";
    cout << myfile;

    cout << "Extracting zf curve attributes from " << myfile << "\n";
    
    string line;
    ifstream zfstats(myfile);
    if (zfstats.is_open()) {
        while ( zfstats.good()) {
            getline(zfstats, line);
            stringstream sstr2(line);
            double val;  
            int j = 0;

            while( sstr2 >> val) {
                    //cout << val << "\n";
                xreal[j++] = val;  /* stored as zmax, z0, q, z10, fmax, fhalfwidth */
            }
        }    
    }
    zfstats.close();        
}

void Individual::copy(Individual *ind1) {
    this->rank = ind1->rank;
    this->constr_violation = ind1->constr_violation;
    this->crowding_distance = ind1->crowding_distance;
    for (int i = 0; i < NUM_PARS; i++) {
        this->vars[i] = ind1->vars[i];
    }
    for (int j = 0; j < NUM_OBJS; j++) {
        this->objs[j] = ind1->objs[j];
    }
    for (int k = 0; k < NUM_CON; k++) {
        this->constr[k] = ind1->constr[k];
    }
    this->ndom = 0;
    this->ss.clear();
}


