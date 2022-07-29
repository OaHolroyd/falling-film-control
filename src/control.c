#include "control.h"

#include "control-internals.h"
#include "control-pair.h"

#include "c-utils.h"


enum CONTROL_TYPE CT;

/* ========================================================================== */
/*   STRATEGY-SPECIFIC FUNCTION DECLARATIONS                                  */
/* ========================================================================== */
/* sets up the specific control strategy */
void (*s_set)();

/* frees varaibles specific to the selected strategy */
void (*s_free)();

/* steps the specific control system forward in time given the interfacial
   height */
void (*control_step)(double *h);

/* returns the estimator as a function of x */
double (*estimator)(double x);


/* ========================================================================== */
/*   FUNCTION DEFINITIONS                                                     */
/* ========================================================================== */
/* returns the baseplate control velocity as a function of x */
double control(double x) {
  double vc = 0.0;

  for (int i = 0; i < M; i++) {
    vc += Amag[i] * actuator(x-Aloc[i]);
  } // i end

  return -ALPHA*vc;
}

/* computes the incremental cost from the current timestep from the interface */
double control_cost(double *h) {
  double cost = 0.0;

  /* control cost */
  // TODO: is this simplification correct or do we need the full integral?
  for (int i = 0; i < M; i++) {
    cost += (1-MU) * Amag[i]*Amag[i];
  } // i end

  /* interfacial cost */
  for (int i = 0; i < N; i++) {
    cost += MU * (h[i]-1)*(h[i]-1);
  } // i end

  return cost;
}

/* set up the control system */
// TODO: explain parameters properly
void control_set(enum CONTROL_TYPE ct, int m, int p, double w, double alpha, double mu, double del, double lx, int n) {
  /* control strategy independent setup */
  internal_control_set(m, p, w, alpha, mu, del, lx, n);

  /* set strategy specific functions */
  CT = ct;
  switch (CT) {
    case PAIR:
      s_set = &pair_set;
      s_free = &pair_free;
      control_step = &pair_step;
      estimator = &pair_estimator;
      break;
    case STATIC:
      ABORT("static control type not implemented yet");
      break;
    case DYNAMIC:
      ABORT("dynamic control type not implemented yet");
      break;
    default :
      ABORT("invalid control type %d", ct);
  }

  /* call strategy specific setup */
  s_set(m, p, w, alpha, mu, del, lx, n);
}

/* frees the control system */
void control_free() {
  internal_control_free();
}
