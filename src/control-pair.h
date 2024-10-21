#ifndef CONTROL_PAIR_H
#define CONTROL_PAIR_H

#include "c-utils.h"
#include "control-core.h"


/* ========================================================================== */
/*   FUNCTION DEFINITIONS                                                     */
/* ========================================================================== */
/* [REQUIRED] internal setup */
void pair_set(void) {
  /* pair requires M and P to be the same */
  if (M != P) {
    ABORT("M and P must be equal for paired controls");
  }
}

/* [REQUIRED] internal free */
void pair_free(void) {}

/* [REQUIRED] steps the system forward in time given the interfacial height */
void pair_step(double dt, double *h, double *q) {
  for (int i = 0; i < M; i++) {
    Amag[i] = ALPHA*(interp(Aloc[i]-DEL, h) - 1.0);
  } // i end
}

/* [REQUIRED] returns the estimator as a function of x */
double pair_estimator(double x) {

  return 0.0;
}

/* [REQUIRED] outputs the internal matrices */
void pair_output(void) {
  // TODO: work out matrix
}

/* [REQUIRED] generates the control matrix CM = a*F*Phi */
void pair_matrix(double **CM) {
  /* forcing matrix */
  double **F = malloc_f2d(N, M);
  forcing_matrix(F);

  /* observer matrix (actually the transpose) */
  double **Phi = malloc_f2d(N, P);
  benney_observer(Phi);

  for (int i = 0; i < N; i++) {
    for (int j = 0; j < N; j++) {
      CM[i][j] = 0.0;
      for (int k = 0; k < M; k++) {
        CM[i][j] += ALPHA*F[i][k]*Phi[j][k];
      } // k end
    } // j end
  } // i end

  free_2d(F);
}


#endif
