#ifndef CONTROL_LQR_H
#define CONTROL_LQR_H

#include <math.h>

#include "c-utils.h"
#include "linalg.h"
#include "control-core.h"


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
  double **A = malloc_f2d(2*N, 2*N);
  wr_jacobian(A);

  /* actuator matrix */
  double **B = malloc_f2d(2*N, M);
  wr_actuator(B);

  /* full control matrix */
  dlqr(A, B, DX*MU, 1-MU, 2*N, M, lqr_k);

  free_2d(A);
  free_2d(B);
}


/* ========================================================================== */
/*   FUNCTION DEFINITIONS                                                     */
/* ========================================================================== */
/* [REQUIRED] internal setup */
void lqr_set(void) {
  /* pick from the available ROMs */
  switch (RT) {
    case BENNEY:
      LQR_K = malloc_f2d(M, N);
      lqr_benney_compute_K(LQR_K);
      break;
    case WR:
      LQR_K = malloc_f2d(M, 2*N);
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
void lqr_step(double dt, double *h, double *q) {
  /* f = K * (h-1) */
  for (int i = 0; i < M; i++) {
    Amag[i] = 0.0;
    for (int j = 0; j < N; j++) {
      Amag[i] += LQR_K[i][j] * (interp(ITOX(j), h) - 1.0);
    } // j end

    for (int j = 0; j < N; j++) {
      Amag[i] += LQR_K[i][N+j] * (interp(ITOX(j), q) - 2.0/3.0);
    } // j end
  } // i end
}

/* [REQUIRED] returns the estimator as a function of x */
double lqr_estimator(double x) {

  return 0.0;
}

/* [REQUIRED] outputs the internal matrices */
void lqr_output(void) {
  if (RT == BENNEY) {
    output_d2d("out/K.dat", LQR_K, M, N);
  } else {
    output_d2d("out/K.dat", LQR_K, M, 2*N);
  }
}

/* [REQUIRED] generates the control matrix CM = F*K */
void lqr_matrix(double **CM) {
  /* forcing matrix */
  double **F = malloc_f2d(N, M);
  forcing_matrix(F);

  for (int i = 0; i < N; i++) {
    for (int j = 0; j < 2 * N; j++) {
      CM[i][j] = 0.0;
      for (int k = 0; k < M; k++) {
        CM[i][j] += F[i][k]*LQR_K[k][j];
      } // k end
    } // j end
  } // i end

  free_2d(F);
}

#endif
