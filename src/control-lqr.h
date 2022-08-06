#ifndef CONTROL_LQR_H
#define CONTROL_LQR_H

#include <math.h>

#include "c-utils.h"
#include "linalg.h"
#include "control-internals.h"


static double **LQR_K; /* control operator */


/* ========================================================================== */
/*   AUXILIARY FUNCTION DEFINITIONS                                           */
/* ========================================================================== */
/* compute the control matrix in the Benney case */
void lqr_benney_compute_K(double **lqr_k) {
  /* Jacobian */
  double **J = malloc_f2d(N, N);
  benney_jacobian(J);

  /* actuator matrix */
  double **Psi = malloc_f2d(N, M);
  benney_actuator(Psi);

  /* control matrix */
  dlqr(J, Psi, DX*MU, 1-MU, N, M, lqr_k);

  free_2d(J);
  free_2d(Psi);
}

/* compute the control matrix in the weighted-residuals case */
void lqr_wr_compute_K(double **lqr_k) {
    /* Jacobian */
  double **J = malloc_f2d(2*N, 2*N);
  wr_jacobian(J);

  /* actuator matrix */
  double **Psi = malloc_f2d(2*N, M);
  wr_actuator(Psi);

  /* full control matrix */
  double **K = malloc_f2d(M, 2*N);
  dlqr(J, Psi, DX*MU, 1-MU, 2*N, M, K);

  /* apply flux approximation */
  for (int i = 0; i < M; i++) {
    for (int j = 0; j < N; j++) {
      lqr_k[i][j] = K[i][j] + 2/3 * K[i][j+N];
    } // j end
  } // i end

  free_2d(J);
  free_2d(Psi);
  free_2d(K);
}


/* ========================================================================== */
/*   FUNCTION DEFINITIONS                                                     */
/* ========================================================================== */
/* [REQUIRED] internal setup */
void lqr_set(void) {
  LQR_K = malloc_f2d(M, N);

  /* pick from the available ROMs */
  switch (RT) {
    case BENNEY:
      lqr_benney_compute_K(LQR_K);
      break;
    case WR:
      lqr_wr_compute_K(LQR_K);
      break;
    default :
      ABORT("invalid ROM type %d", RT);
  }
}

/* [REQUIRED] internal free */
void lqr_free(void) {

  free(LQR_K);
}

/* [REQUIRED] steps the system forward in time given the interfacial height */
void lqr_step(double dt, double *h) {
  /* f = K * (h-1) */
  for (int i = 0; i < M; i++) {
    Amag[i] = 0.0;
    for (int j = 0; j < N; j++) {
      Amag[i] += LQR_K[i][j] * (interp(ITOX(j), h) - 1.0);
    } // j end
  } // i end
}

/* [REQUIRED] returns the estimator as a function of x */
double lqr_estimator(double x) {

  return 0.0;
}


#endif
