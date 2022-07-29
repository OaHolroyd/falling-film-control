#ifndef CONTROL_PAIR_H
#define CONTROL_PAIR_H

#include "c-utils.h"
#include "control-internals.h"


/* ========================================================================== */
/*   FUNCTION DEFINITIONS                                                     */
/* ========================================================================== */
/* [REQUIRED] internal setup */
void pair_set() {
  /* pair requires M and P to be the same */
  if (M != P) {
    ABORT("M and P must be equal for paired controls");
  }
}

/* [REQUIRED] internal free */
void pair_free() {}

/* [REQUIRED] steps the system forward in time given the interfacial height */
void pair_step(double *h) {
  for (int i = 0; i < M; i++) {
    Amag[i] = interp(Aloc[i]-DEL, h) - 1.0;
  } // i end
}

/* [REQUIRED] returns the estimator as a function of x */
double pair_estimator(double x) {

  return 0.0;
}

#endif
