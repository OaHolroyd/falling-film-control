#include "control.h"

#include "control-internals.h"
#include "control-pair.h"
#include "control-static.h"
#include "control-dynamic.h"

#include "c-utils.h"


static control_t CT;


/* ========================================================================== */
/*   STRATEGY-SPECIFIC FUNCTION DECLARATIONS                                  */
/* ========================================================================== */
/* sets up the specific control strategy */
void (*s_set)(void);

/* frees varaibles specific to the selected strategy */
void (*s_free)(void);

/* steps the specific control system forward in time given the interfacial
   height */
void (*control_step)(double dt, double *h);

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

  double xi, ai;
  for (int i = 0; i < N; i++) {
    xi = ITOX(i);

    /* control cost */
    ai = actuator(xi);
    cost += (1-MU) * ai*ai;

    /* interfacial cost */
    cost += MU * (h[i]-1)*(h[i]-1);
  } // i end

  /* integral scaling */
  cost *= DX;

  return cost;
}

/* set up the control system */
// TODO: explain parameters properly
void control_set(control_t ct, rom_t rt, int m, int p, double w, double alpha, double mu, double del, double lx, int n, double re, double ca, double theta) {
  /* control strategy independent setup */
  internal_control_set(rt, m, p, w, alpha, mu, del, lx, n, re, ca, theta);

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
      s_set = &static_set;
      s_free = &static_free;
      control_step = &static_step;
      estimator = &static_estimator;
      break;
    case DYNAMIC:
      s_set = &dynamic_set;
      s_free = &dynamic_free;
      control_step = &dynamic_step;
      estimator = &dynamic_estimator;
      break;
    default :
      ABORT("invalid control type %d", ct);
  }

  /* call strategy specific setup */
  s_set();
}

/* frees the control system */
void control_free(void) {
  internal_control_free();
  s_free();
}
