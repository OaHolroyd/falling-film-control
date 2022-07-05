#ifndef CONTROL_PAIR_H
#define CONTROL_PAIR_H

#include <math.h>

#include "parallel.h"
#include "film-utils.h"
#include "params.h"
#include "control.h"


/* ========================================================================== */
/*   FUNCTION DEFINITIONS                                                     */
/* ========================================================================== */
/* [REQUIRED] Sets the control variables (call after set_params) */
void set_Cparams() {

  internal_set_Cparams();
}

/* [REQUIRED] frees control variables */
void control_free() {

  internal_control_free();
}

/* [REQUIRED] sets the control magnitudes and returns the cost */
void control_set_magnitudes() {
  /* basic paired observer-actuators */
  for (int i = 0; i < C_M; i++) {
    C_mag[i] = interfacial_height(C_loc[i] - C_PHI) - 1;
  } // i end
}


#endif
