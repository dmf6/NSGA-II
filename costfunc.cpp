#include <cstdlib>
#include <math.h>
#include <iostream>
#include "globals.h"

void costfunc(double *x, double *obj, double *constr) {
    obj[0] = pow(x[0], 2);
    obj[1] = pow(x[0] - 2, 2);

    /*ZDT4*/
    // double f1, f2, g, h;
    // int i;
    // f1 = x[0];
    // g = 0.0;
    // for (i=1; i<10; i++)
    // {
    //     g += x[i]*x[i] - 10.0*cos(4.0*PI*x[i]);
    // }
    // g += 91.0;
    // h = 1.0 - sqrt(f1/g);
    // f2 = g*h;
    // obj[0] = f1;
    // obj[1] = f2;

    
    // double f1, f2, g, h;
  
    // f1 = 1.0 - ( exp(-4.0*x[0]) ) * pow( (sin(6.0*PI*x[0])),6.0 );
    // g = 0.0;
    // for (int i=1; i<10; i++)
    // {
    //     g += x[i];
    // }
    // g = g/9.0;
    // g = pow(g,0.25);
    // g = 1.0 + 9.0*g;
    // h = 1.0 - pow((f1/g),2.0);
    // f2 = g*h;
    // obj[0] = f1;
    // obj[1] = f2;
}
