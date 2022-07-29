#ifndef CONTROL_INTERNALS_H
#define CONTROL_INTERNALS_H

#include <stdlib.h>
#include <math.h>
#include <complex.h>

#include "control.h"


/* convert from index to location */
#define ITOX(i) (DX*(i+0.5))

/* convert from location to index */
#define XTOI(x) ((int)(x/DX - 0.5))


/* ========================================================================== */
/*   GLOBAL DEFINITIONS                                                       */
/* ========================================================================== */
/* film constants */
static int N; // number of gridcells
static double LX; // domain length
static double DX; // gridspacing
static double RE; // Reynolds' number
static double CA; // capillary number
static double THETA; // plate angle

/* control constants */
static int M; // number of actuators
static int P; // number of observers
static double W; // width parameter
static double ALPHA; // control strength parameter
static double MU; // control cost parameter
static double NORM; // normalising constant
static double DEL; // observer offset (upstream)
static rom_t RT; // type of reduced order model

/* location arrays */
static double *Aloc; // actuator locations
static double *Oloc; // observer locations

/* variable arrays */
static double *Amag; // actuator magnitudes


/* ========================================================================== */
/*   FUNCTION DEFINITIONS                                                     */
/* ========================================================================== */
/* interpolates the discretised film height h at a position x */
double interp(double x, double *h) {
  /* account for periodicity */
  if (x < 0) {
    x += LX;
  } else if (x > LX) {
    x -= LX;
  }

  /* get neighbouring indices */
  double i = x/DX - 0.5;
  int i0 = floor(i);
  int i1 = ceil(i);

  /* cover the case that x is on a gridpoint */
  if (i0 == i1) {
    return h[i0];
  }

  /* linear interpolation */
  return h[i0]*(i1-i) + h[i1]*(i-i0);
}

/* (periodic) actuator function (only valid on [-3Lx/2, 3Lx/2]) */
double actuator(double x) {
  /* account for periodicity */
  if (x < -LX/2) {
    x += LX;
  } else if (x > LX/2) {
    x -= LX;
  }

  return NORM*exp((cos(2*M_PI*x/LX)-1.0)/(W*W));
}

/* sets the common control parameters and allocates common memory */
void internal_control_set(rom_t rt, int m, int p, double w, double alpha, double mu, double del, double lx, int n, double re, double ca, double theta) {
  /* set constants */
  N = n;
  LX = lx;
  DX = LX/N;
  RE = re;
  CA = ca;
  THETA = theta;

  M = m;
  P = p;
  W = w;
  ALPHA = alpha;
  MU = mu;
  DEL = del;
  RT = rt;

  /* actuator locations/magnitudes */
  Aloc = malloc(M*sizeof(double));
  Amag = malloc(M*sizeof(double));
  for (int i = 0; i < M; i++) {
    Aloc[i] = (i+0.5) * LX / M;
    Amag[i] = 0.0;
  } // i end

  /* observer locations */
  Oloc = malloc(P*sizeof(double));
  for (int i = 0; i < P; i++) {
    Oloc[i] = (i+0.5) * LX / P - DEL;

    /* wrap via periodicity */
    // TODO: this will break if |DEL| > LX
    if (Oloc[i] < 0) {
      Oloc[i] += LX;
    } else if (Oloc[i] > LX) {
      Oloc[i] -= LX;
    }
  } // i end

  /* control normaliser */
  NORM = 1.0;
  double integral = 0.0;
  for (int i = 0; i < N; i++) {
    integral += actuator(DX*i - LX/2);
  } // i end
  NORM = 1.0/(DX*integral);
}

/* frees the common memory */
void internal_control_free() {
  free(Aloc);
  free(Oloc);
  free(Amag);
}


#endif
