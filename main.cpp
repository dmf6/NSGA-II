#include <iostream>
#include <cstdlib>
#include <stdio.h>
#include <fstream>
#include <sstream>
#include <cstring>
#include <qapplication.h>
#include <QtGui>
#include <qgraphicsitem.h>
#include "population.h"
#include "window.h"

double* Population::min_var = new double[NUM_PARS];
double* Population::max_var = new double[NUM_PARS];

int main(int argc, char *argv[]) {
        /* GUI STUFF */
    QApplication app( argc, argv );
    
    Window widget;
    widget.setWindowTitle(QT_TRANSLATE_NOOP(QGraphicsView, "My Qt drawing canvas"));
    widget.drawRectangle();
    
    widget.show();
         
     
    Population *parent_pop = new Population(POP_SIZE, NUM_PARS, NUM_OBJS, NUM_CON);
    Population *child_pop = new Population(POP_SIZE, NUM_PARS, NUM_OBJS, NUM_CON);
    Population *mixed_pop = new Population(2*POP_SIZE, NUM_PARS, NUM_OBJS, NUM_CON);
    ifstream limitsfile("limits.dat");
    parent_pop->readLimits(limitsfile);
    limitsfile.close();
    
    //     /* initialize parent population with randomly generated
    //      * numbers using boost libraries */
    parent_pop->initialize();
    //     /* send initial population to a file */
    ofstream os;
    os.open("initial_pop.dat");
    parent_pop->print_pop(os);
    os.close();

        /* initial parallel simulations */
        //parent_pop->simulate();
        // evaluate entire population
    parent_pop->evaluate();
    // //     /* perform non-dominated sorting and assign crowding distance */
    
    parent_pop->nondominated_sort(); 

    //     /*  Start evolution */
      for (int i = 1; i < NUM_GEN; i++) {
              // cout << endl << "Starting generation " << i << "\n";
     /* select parents and produce offspring by SBX crossover */
          parent_pop->select(child_pop);
       /*mutate children */
          child_pop->mutate_pop();
            /* this is done in parallel */
             //child_pop->simulate();
             /* evaluate how good of a match to target zf profile */
          child_pop->evaluate();
           /* populate mixed_pop with parent and child populations */
          mixed_pop->merge(parent_pop, child_pop);
           /*
            *  Sort population into different nondomination levels. Each
            *  individual is assigned a fitness equal to its
            *  non-domination level (1 being the best)
            */
          mixed_pop->nondominated_sort();  //non-dominated sort on all fronts
          //  /* create new parent population from mixed population */
          mixed_pop->generateNewParent(parent_pop);
      }

          //cout << "Generations finished, now sending solutions to final_pop.dat" << endl;


    /* print parameters */
    ofstream os1;
    os1.open("final_pop.dat");
    parent_pop->print_pop(os1);
    os1.close();
        /* print objectives */
     ofstream os2;
    os2.open("final_pop_objectives.dat");
    parent_pop->printObjectives(os2);
    os2.close();
    
    // delete parent_pop;
    // delete child_pop;
    // delete mixed_pop;

    return app.exec();
}
