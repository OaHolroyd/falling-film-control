#ifndef CONTROL_H
#define CONTROL_H

#include <math.h>

#include "film-utils.h"
#include "params.h"

double C_cost = 0.0; // total control cost

/* ========================================================================== */
/*   FUNCTION DEFINITIONS                                                     */
/* ========================================================================== */
/* discrete actuator function */
double actuator(double x) {

  return C_norm*exp((cos(2*M_PI*x/LX)-1.0)/(C_W*C_W));
}

/* Sets the control variables (call after set_params) */
void set_Cparams() {
  C_loc = malloc(C_M*sizeof(double));
  C_mag = malloc(C_M*sizeof(double));

  /* control locations */
  double dc = LX/C_M;
  for (int i = 0; i < C_M; i++) {
    C_loc[i] = (i+0.5)*dc;
    C_mag[i] = 0.0;
  } // i end

  /* set control normaliser */
  C_norm = 1.0;
  double I = 0.0;
  int nx = (1<<LEVEL);
  double dx = LX/((double)(nx));
  for (int i = 0; i < nx; i++) {
    I += actuator(dx*i - LX/2.0);
  } // i end
  C_norm = 1.0/(dx*I);
}

/* returns the baseplate control velocity as a function of x */
double control(double x) {
  double vc = 0.0;

  for (int i = 0; i < C_M; i++) {
    vc += C_mag[i] * actuator(x-C_loc[i]);
  } // i end

  return -C_ALPHA*vc;
}

/* sets the control magnitudes and returns the cost */
double control_set_magnitudes() {
  double c = 0.0;

  for (int i = 0; i < C_M; i++) {
    C_mag[i] = interfacial_height(C_loc[i] - C_PHI) - 1;
    c += C_mag[i]*C_mag[i];
  } // i end

  return c;
}


#endif
