#ifndef CONTROL_QQR_H
#define CONTROL_QQR_H

#include <math.h>

#include "c-utils.h"
#include "linalg.h"
#include "control-core.h"


static double **QQR_C0; /* control operator */
static double **QQR_C1; /* control operator */


/* ========================================================================== */
/*   AUXILIARY FUNCTION DEFINITIONS                                           */
/* ========================================================================== */
/* compute the control matrix in the Benney case */
void qqr_benney_compute_matrices(double **qqr_c0, double **qqr_c1) {
  /* Jacobian */
  double **J = malloc_f2d(N, N);
  benney_jacobian(J);

  /* actuator matrix */
  double **Psi = malloc_f2d(N, M);
  benney_actuator(Psi);

  /* control matrix */
  dlqr(J, Psi, DX*MU, 1-MU, N, M, qqr_c0);

  free_2d(J);
  free_2d(Psi);
}

/* compute the control matrix in the weighted-residuals case */
void qqr_wr_compute_matrices(double **qqr_c0, double **qqr_c1) {
  /* Jacobian */
  double **A = malloc_f2d(2*N, 2*N);
  wr_jacobian(A);

  /* actuator matrix */
  double **B = malloc_f2d(2*N, M);
  wr_actuator(B);

  /* full control matrix */
  dlqr(A, B, DX*MU, 1-MU, 2*N, M, qqr_c0);

  free_2d(A);
  free_2d(B);
}


/* ========================================================================== */
/*   FUNCTION DEFINITIONS                                                     */
/* ========================================================================== */
/* [REQUIRED] internal setup */
void qqr_set(void) {
  QQR_C0 = malloc_f2d(M, 2*N);
  QQR_C1 = malloc_f2d(M, 2*N);

  /* pick from the available ROMs */
  switch (RT) {
    case BENNEY:
      qqr_benney_compute_matrices(QQR_C0, QQR_C1);
      break;
    case WR:
      qqr_wr_compute_matrices(QQR_C0, QQR_C1);
      break;
    default :
      ABORT("invalid ROM type %d", RT);
  }
}

/* [REQUIRED] internal free */
void qqr_free(void) {
  free(QQR_C0);
  free(QQR_C1);
}

/* [REQUIRED] steps the system forward in time given the interfacial height */
void qqr_step(double dt, double *h, double *q) {
  /* f = K * (h-1) */
  for (int i = 0; i < M; i++) {
    Amag[i] = 0.0;

    /* height component */
    for (int j = 0; j < N; j++) {
      // TODO: compute nonlinear bit correctly
      const double hj = interp(ITOX(j), h) - 1.0;
      const double hhj = 0.0;

      Amag[i] += QQR_C0[i][j] * hj + QQR_C1[i][j] * hhj;
    } // j end

    /* flux component */
    for (int j = 0; j < N; j++) {
      // TODO: compute nonlinear bit correctly
      const double qj = interp(ITOX(j), q) - 2.0/3.0;
      const double qqj = 0.0;

      Amag[i] += QQR_C0[i][N+j] * qj + QQR_C1[i][j] * qqj;
    } // j end
  } // i end
}

/* [REQUIRED] returns the estimator as a function of x */
double qqr_estimator(double x) {

  return 0.0;
}

/* [REQUIRED] outputs the internal matrices */
void qqr_output(void) {
  output_d2d("out/C0.dat", QQR_C0, M, 2*N);
  output_d2d("out/C1.dat", QQR_C1, M, 2*N);
}

/* [REQUIRED] generates the control matrix CM = F*K */
void qqr_matrix(double **CM) {
  exit(1);

  /* forcing matrix */
  double **F = malloc_f2d(N, M);
  forcing_matrix(F);

  for (int i = 0; i < N; i++) {
    for (int j = 0; j < N; j++) {
      CM[i][j] = 0.0;
      for (int k = 0; k < M; k++) {
        CM[i][j] += F[i][k]*QQR_C0[k][j];
      } // k end
    } // j end
  } // i end

  free_2d(F);
}

#endif
