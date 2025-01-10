#ifndef CONTROL_SDRE_H
#define CONTROL_SDRE_H

#include <math.h>

#include "c-utils.h"
#include "linalg.h"
#include "control-core.h"


static double **SDRE_K1; // first-order gain operator
static double **SDRE_K2; // basis for second-order gain correction


/* ========================================================================== */
/*   AUXILIARY FUNCTION DEFINITIONS                                           */
/* ========================================================================== */
/* compute the control matrix in the Benney case */
void sdre_benney_compute_K(void) {
  /* Jacobian */
  double **J = malloc_f2d(N, N);
  benney_jacobian(J);

  /* actuator matrix */
  double **Psi = malloc_f2d(N, M);
  benney_actuator(Psi);

  /* control matrix */
  dlqr(J, Psi, DX*MU, 1-MU, N, M, SDRE_K1);

  free_2d(J);
  free_2d(Psi);
}

/* compute the control matrix in the weighted-residuals case */
void sdre_wr_compute_K(void) {
    /* Jacobian */
  double **A = malloc_f2d(2*N, 2*N);
  wr_jacobian(A);

  /* actuator matrix */
  double **B = malloc_f2d(2*N, M);
  wr_actuator(B);

  /* full control matrix */
  dlqr(A, B, DX*MU, 1-MU, 2*N, M, SDRE_K1);

  free_2d(A);
  free_2d(B);
}


/* ========================================================================== */
/*   FUNCTION DEFINITIONS                                                     */
/* ========================================================================== */
/* [REQUIRED] internal setup */
void sdre_set(void) {
  /* pick from the available ROMs */
  switch (RT) {
    case BENNEY:
      SDRE_K1 = malloc_f2d(M, N);
      sdre_benney_compute_K();
      break;
    case WR:
      fprintf(stderr, "\nCan't do SDRE with WR\n");
      exit(EXIT_SUCCESS);
      SDRE_K1 = malloc_f2d(M, 2*N);
      sdre_wr_compute_K();
      break;
    default :
      ABORT("invalid ROM type %d", RT);
  }
}

/* [REQUIRED] internal free */
void sdre_free(void) {

  free(SDRE_K1);
}

/* [REQUIRED] steps the system forward in time given the interfacial height */
void sdre_step(double dt, double *h, double *q) {
  // RECOMPUTE K, accounting for nonlinearities

  /* f = K * (h-1) */
  for (int i = 0; i < M; i++) {
    Amag[i] = 0.0;
    for (int j = 0; j < N; j++) {
      Amag[i] += SDRE_K1[i][j] * (interp(ITOX(j), h) - 1.0);
    } // j end

    /* only WR uses the flux */
    if (RT == WR) {
      for (int j = 0; j < N; j++) {
        Amag[i] += SDRE_K1[i][N+j] * (interp(ITOX(j), q) - 2.0/3.0);
      } // j end
    }
  } // i end
}

/* [REQUIRED] returns the estimator as a function of x */
double sdre_estimator(double x) {

  return 0.0;
}

/* [REQUIRED] outputs the internal matrices */
void sdre_output(void) {
  if (RT == BENNEY) {
    output_d2d("out/K.dat", SDRE_K1, M, N);
  } else {
    output_d2d("out/K.dat", SDRE_K1, M, 2*N);
  }
}

/* [REQUIRED] generates the control matrix CM = F*K */
void sdre_matrix(double **CM) {
  /* forcing matrix */
  double **F = malloc_f2d(N, M);
  forcing_matrix(F);

  for (int i = 0; i < N; i++) {
    for (int j = 0; j < 2 * N; j++) {
      CM[i][j] = 0.0;
      for (int k = 0; k < M; k++) {
        CM[i][j] += F[i][k]*SDRE_K1[k][j];
      } // k end
    } // j end
  } // i end

  free_2d(F);
}

#endif
