#ifndef CONTROL_H
#define CONTROL_H

#include <math.h>

#include "parallel.h"
#include "film-utils.h"
#include "params.h"


double C_cost = 0.0; // total control cost
double C_norm; // actuator normaliser
double *C_loc; // actuator locations
double *C_mag; // current control magnitudes

/* ========================================================================== */
/*   FUNCTION DEFINITIONS                                                     */
/* ========================================================================== */
/* (periodic) actuator function */
double actuator(double x) {
  /* account for periodicity */
  if (x < -LX/2) {
    x += LX;
  } else if (x > LX/2) {
    x -= LX;
  }

  return C_norm*exp((cos(2*M_PI*x/LX)-1.0)/(C_W*C_W));
}

/* returns the baseplate control velocity as a function of x */
double control(double x) {
  double vc = 0.0;

  for (int i = 0; i < C_M; i++) {
    vc += C_mag[i] * actuator(x-C_loc[i]);
  } // i end

  return -C_ALPHA*vc;
}

/* computes the incremental cost at the current timestep */
double control_cost() {
  /* control cost */
  for (int i = 0; i < C_M; i++) {
    C_cost += dt * (1-C_MU) * C_mag[i] * C_mag[i];
  } // i end

  /* interfacial cost */
  double dh;
  double xi;
  for (int i = 0; i < N; i++) {
    xi = DX*(i+0.5);
    dh = interfacial_height(xi) - 1;
    C_cost += dt * C_MU * dh * dh;
  } // i end
}

/* sets the control variables (call after set_params) */
void internal_set_Cparams() {
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
  double integral = 0.0;
  for (int i = 0; i < N; i++) {
    integral += actuator(DX*i - LX/2.0);
  } // i end
  C_norm = 1.0/(DX*integral);
}

/* frees control variables */
void internal_control_free() {
  free(C_loc);
  free(C_mag);
}


#endif
