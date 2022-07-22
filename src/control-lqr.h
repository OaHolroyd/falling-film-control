#ifndef CONTROL_LQR_H
#define CONTROL_LQR_H

#include <math.h>

#include "c-utils.h"
#include "parallel.h"
#include "film-utils.h"
#include "params.h"
#include "control.h"
#include "lqr.h"


double **C_K; /* control operator */


/* ========================================================================== */
/*   AUXILIARY FUNCTION DEFINITIONS                                           */
/* ========================================================================== */
/* set the Jacobian matrix for the Benney system */
void fill_jacobian_benney(double **J) {
  /* ensure filled with zeros */
  for (int i = 0; i < N; i++) {
    for (int j = 0; j < N; j++) {
      J[i][j] = 0.0;
    } // j end
  } // i end

  double c0 = -1/(3*CA) * (1/(DX*DX*DX*DX));
  double c1 = 1/DX + (2/tan(THETA)/3 - 8*RE/15) * (1/(DX*DX)) + 1/(3*CA) * (4/(DX*DX*DX*DX));
  double c2 = (2/tan(THETA)/3 - 8*RE/15) * (-2/(DX*DX)) - 1/(3*CA) * (6/(DX*DX*DX*DX));
  double c3 = -1/DX + (2/tan(THETA)/3 - 8*RE/15) * (1/(DX*DX)) + 1/(3*CA) * (4/(DX*DX*DX*DX));
  double c4 = -1/(3*CA) * (1/(DX*DX*DX*DX));
  for (int i = 2; i < N-2; i++) {
    J[i][i-2] = c0;
    J[i][i-1] = c1;
    J[i][i] = c2;
    J[i][i+1] = c3;
    J[i][i+2] = c4;
  } // i end

  J[0][N-2] = c0;
  J[0][N-1] = c1;
  J[0][0] = c2;
  J[0][1] = c3;
  J[0][2] = c4;

  J[1][N-1] = c0;
  J[1][0] = c1;
  J[1][1] = c2;
  J[1][1+1] = c3;
  J[1][1+2] = c4;

  J[N-2][N-4] = c0;
  J[N-2][N-3] = c1;
  J[N-2][N-2] = c2;
  J[N-2][N-1] = c3;
  J[N-2][0] = c4;

  J[N-1][N-3] = c0;
  J[N-1][N-2] = c1;
  J[N-1][N-1] = c2;
  J[N-1][0] = c3;
  J[N-1][1] = c4;
}

/* set the discrete actuator matrix for the Benney system */
void fill_actuators_benney(double **Psi) {
  /* forcing matrix (TODO: see if rotating this makes it faster) */
  double **F = malloc_f2d(N, C_M);
  for (int i = 0; i < N; i++) {
    for (int j = 0; j < C_M; j++) {
      F[i][j] = actuator(DX*(i+0.5)-C_loc[j]);
    } // j end
  } // i end

  /* actuator matrix (which is equal to the forcing matrix in the simplified
     Benney case) */
  for (int i = 1; i < N-1; i++) {
    for (int j = 0; j < C_M; j++) {
      Psi[i][j] = F[i][j] + (RE/(3*DX)) * (F[i+1][j]-F[i-1][j]);
    } // j end
  } // i end
  for (int j = 0; j < C_M; j++) {
    Psi[0][j] = F[0][j] + (RE/(3*DX)) * (F[1][j]-F[N-1][j]);
    Psi[N-1][j] = F[N-1][j] + (RE/(3*DX)) * (F[0][j]-F[N-2][j]);
  } // j end

  free_2d((void **)F);
}

/* compute the control matrix for the Benney system */
void compute_K_benney() {
  C_K = malloc_f2d(C_M, N);

  /* interim arrays */
  double **J = malloc_f2d(N, N);
  fill_jacobian_benney(J);

  double **Psi = malloc_f2d(N, C_M);
  fill_actuators_benney(Psi);

  /* compute K */
  lqr(J, Psi, C_MU, 1-C_MU, N, C_M, C_K);

  free_2d((void **)J);
  free_2d((void **)Psi);
}

/* set the Jacobian matrix for the WR system */
void fill_jacobian_wr(double **J) {
  /* ensure filled with zeros */
  for (int i = 0; i < 2*N; i++) {
    for (int j = 0; j < 2*N; j++) {
      J[i][j] = 0.0;
    } // j end
  } // i end

  double c0, c1, c2, c3, c4;

  /* top right */
  c1 = 1/(2*DX);
  c3 = -1/(2*DX);
  for (int i = 1; i < N-1; i++) {
    J[i][N+i-1] = c1;
    J[i][N+i+1] = c3;
  } // i end

  J[0][N+N-1] = c1;
  J[0][N+1] = c3;

  J[N-1][N+N-2] = c1;
  J[N-1][N+0] = c3;


  /* bottom left */
  c0 = 1/(3*CA) * (-1/(2*DX*DX*DX));
  c1 = (8*RE/35 - 2/tan(THETA)/3) * (-1/(2*DX)) + 1/(3*CA) * (1/(DX*DX*DX));
  c2 = 2;
  c3 = (8*RE/35 - 2/tan(THETA)/3) * (1/(2*DX)) + 1/(3*CA) * (-1/(DX*DX*DX));
  c4 = 1/(3*CA) * (1/(2*DX*DX*DX));
  for (int i = 2; i < N-2; i++) {
    J[N+i][i-2] = c0;
    J[N+i][i-1] = c1;
    J[N+i][i] = c2;
    J[N+i][i+1] = c3;
    J[N+i][i+2] = c4;
  } // i end

  J[N+0][N-2] = c0;
  J[N+0][N-1] = c1;
  J[N+0][0] = c2;
  J[N+0][1] = c3;
  J[N+0][2] = c4;

  J[N+1][N-1] = c0;
  J[N+1][0] = c1;
  J[N+1][1] = c2;
  J[N+1][1+1] = c3;
  J[N+1][1+2] = c4;

  J[N+N-2][N-4] = c0;
  J[N+N-2][N-3] = c1;
  J[N+N-2][N-2] = c2;
  J[N+N-2][N-1] = c3;
  J[N+N-2][0] = c4;

  J[N+N-1][N-3] = c0;
  J[N+N-1][N-2] = c1;
  J[N+N-1][N-1] = c2;
  J[N+N-1][0] = c3;
  J[N+N-1][1] = c4;


  /* bottom right */
  c1 = 6*RE/105 * 1/(2*DX);
  c2 = -1;
  c3 = 6*RE/105 * -1/(2*DX);
  for (int i = 1; i < N-1; i++) {
    J[N+i][N+i-1] = c1;
    J[N+i][N+i] = c2;
    J[N+i][N+i+1] = c3;
  } // i end

  J[N+0][N+N-1] = c1;
  J[N+0][N+0] = c2;
  J[N+0][N+1] = c3;

  J[N+N-1][N+N-2] = c1;
  J[N+N-1][N+N-1] = c2;
  J[N+N-1][N+0] = c3;
}

/* set the discrete actuator matrix for the WR system */
void fill_actuators_wr(double **Psi) {
  /* forcing matrix (TODO: see if rotating this makes it faster) */
  double **F = malloc_f2d(N, C_M);
  for (int i = 0; i < N; i++) {
    for (int j = 0; j < C_M; j++) {
      F[i][j] = actuator(DX*(i+0.5)-C_loc[j]);
    } // j end
  } // i end

  /* actuator matrix */
  for (int i = 0; i < N; i++) {
    for (int j = 0; j < C_M; j++) {
      Psi[i][j] = F[i][j];
      Psi[i+N][j] = 2*RE/15 * F[i][j];
    } // j end
  } // i end

  free_2d((void **)F);
}

/* compute the control matrix for the WR system */
void compute_K_wr() {
  double **K = malloc_f2d(C_M, 2*N);

  /* interim arrays */
  double **J = malloc_f2d(2*N, 2*N);
  fill_jacobian_wr(J);

  double **Psi = malloc_f2d(2*N, C_M);
  fill_actuators_wr(Psi);

  /* compute K */
  lqr(J, Psi, C_MU, 1-C_MU, 2*N, C_M, K);

  /* use q = 2/3 h to remove flux observations */
  C_K = malloc_f2d(C_M, N);
  for (int i = 0; i < C_M; i++) {
    for (int j = 0; j < N; j++) {
      C_K[i][j] = K[i][j] + 2/3 * K[i][j+N];
    } // j end
  } // i end

  free_2d((void **)K);
  free_2d((void **)J);
  free_2d((void **)Psi);
}


/* ========================================================================== */
/*   FUNCTION DEFINITIONS                                                     */
/* ========================================================================== */
/* [REQUIRED] sets the control variables (call after set_params) */
void set_Cparams() {
  internal_set_Cparams();

  /* set control operator */
  // TODO: fix issues with K_wr
  compute_K_benney();
  // compute_K_wr();
}

/* [REQUIRED] frees control variables */
void control_free() {
  internal_control_free();

  free_2d((void **)C_K);
}

/* [REQUIRED] sets the control magnitudes and returns the cost */
void control_set_magnitudes() {
  /* basic paired observer-actuators */
  for (int i = 0; i < C_M; i++) {
    C_mag[i] = 0.0;
    for (int j = 0; j < N; j++) {
      C_mag[i] -= C_K[i][j] * (interfacial_height(DX*(j+0.5)) - 1);
    } // j end
  } // i end
}

/* [REQUIRED] returns the film height estimator */
double estimator(double x) {
  /* account for periodicity */
  if (x < 0) {
    x += LX;
  } else if (x > LX) {
    x -= LX;
  }

  /* sum the terms from the Fourier series */
  // TODO: implement
  return 0.0;
}


#endif
