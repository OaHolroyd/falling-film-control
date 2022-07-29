#ifndef CONTROL_H
#define CONTROL_H

#include <complex.h>


enum CONTROL_TYPE { PAIR, STATIC, DYNAMIC };


/* ========================================================================== */
/*   FUNCTION DECLARATIONS                                                    */
/* ========================================================================== */
/* returns the baseplate control velocity as a function of x */
double control(double x);

/* computes the incremental cost from the current timestep from the interface */
double control_cost(double *h);

/* set up the control system */
void control_set(enum CONTROL_TYPE ct, int m, int p, double w, double alpha, double mu, double del, double lx, int n);

/* frees the control system */
void control_free();

/* steps the specific control system forward in time given the interfacial
   height */
extern void (*control_step)(double *h);

/* returns the estimator as a function of x */
extern double (*estimator)(double x);


#endif
