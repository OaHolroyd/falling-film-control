#ifndef CONTROL_LQR_H
#define CONTROL_LQR_H

#include <math.h>

#include "array.h"
#include "parallel.h"
#include "film-utils.h"
#include "params.h"
#include "control.h"
#include "lqr.h"


double **C_K; /* control operator */


/* ========================================================================== */
/*   AUXILIARY FUNCTION DEFINITIONS                                           */
/* ========================================================================== */
/* set the Jacobian matrix */
void fill_jacobian(double **J) {
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

void fill_actuators(double **Psi) {
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


/* ========================================================================== */
/*   FUNCTION DEFINITIONS                                                     */
/* ========================================================================== */
/* [REQUIRED] sets the control variables (call after set_params) */
void set_Cparams() {
  internal_set_Cparams();

  /* control operator (TODO: BENNEY SYSTEM) */
  C_K = malloc_f2d(C_M, N);

  /* interim arrays */
  double **J = malloc_f2d(N, N);
  fill_jacobian(J);

  double **Psi = malloc_f2d(N, C_M);
  fill_actuators(Psi);

  /* compute K */
  lqr(J, Psi, C_MU, 1-C_MU, N, C_M, C_K);

  free_2d((void **)J);
  free_2d((void **)Psi);
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


#endif
