#ifndef CONTROL_ESTIMATOR_H
#define CONTROL_ESTIMATOR_H

#include <math.h>

#include "c-utils.h"
#include "control-core.h"
#include "linalg.h"

static double *EST_h; // height estimate

/* ========================================================================== */
/*   AUXILIARY FUNCTION DEFINITIONS                                           */
/* ========================================================================== */
/**/
void est_update(double dt, double *h) {
  for (int i = 0; i < N; i++) {
    EST_h[i] = 1.0 + 0.8 * (h[i] - 1.0);
  }
}


/* ========================================================================== */
/*   FUNCTION DEFINITIONS                                                     */
/* ========================================================================== */
/* [REQUIRED] internal setup */
void est_set(void) {
  EST_h = malloc(N * sizeof(double));

  /* pick from the available ROMs */
  switch (RT) {
  case BENNEY:
    //      lqr_benney_compute_K(LQR_K);
    break;
  case WR:
    //      lqr_wr_compute_K(LQR_K);
    break;
  default:
    ABORT("invalid ROM type %d", RT);
  }
}

/* [REQUIRED] internal free */
void est_free(void) { free(EST_h); }

/* [REQUIRED] steps the system forward in time given the interfacial height */
void est_step(double dt, double *h) {
  /* update the estimator */
  est_update(dt, h);

  /* f = K * (h-1) */
  // no control
  for (int i = 0; i < M; i++) {
    Amag[i] = 0.0;
    //    for (int j = 0; j < N; j++) {
    //      Amag[i] += LQR_K[i][j] * (interp(ITOX(j), h) - 1.0);
    //    } // j end
  } // i end
}

/* [REQUIRED] returns the estimator as a function of x */
double est_estimator(double x) {
  /* interpolate from EST_h to x */
  return interp(x, EST_h);
}

/* [REQUIRED] outputs the internal matrices */
void est_output(void) {
  //  output_d2d("out/K.dat", LQR_K, M, N);
}

/* [REQUIRED] generates the control matrix CM = F*K */
void est_matrix(double **CM) {
  /* forcing matrix */
  double **F = malloc_f2d(N, M);
  forcing_matrix(F);

  for (int i = 0; i < N; i++) {
    for (int j = 0; j < N; j++) {
      CM[i][j] = 0.0;
      for (int k = 0; k < M; k++) {
        //        CM[i][j] += F[i][k]*LQR_K[k][j];
      } // k end
    } // j end
  } // i end

  free_2d(F);
}

#endif
