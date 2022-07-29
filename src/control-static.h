#ifndef CONTROL_STATIC_H
#define CONTROL_STATIC_H

#include <math.h>

#include "c-utils.h"
#include "lqr.h"
#include "control-internals.h"


static double **STATIC_K; /* control operator */


void benney_compute_K() {
  /* Jacobian */
  double **J = malloc_f2d(N, N);

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


  /* forcing matrix (TODO: see if rotating this makes it faster) */
  double **F = malloc_f2d(N, M);
  for (int i = 0; i < N; i++) {
    for (int j = 0; j < M; j++) {
      F[i][j] = actuator(ITOX(i)-Aloc[j]);
    } // j end
  } // i end


  /* actuator matrix */
  double **Psi = malloc_f2d(N, M);
  for (int i = 1; i < N-1; i++) {
    for (int j = 0; j < M; j++) {
      Psi[i][j] = F[i][j] + (RE/(3*DX)) * (F[i+1][j]-F[i-1][j]);
    } // j end
  } // i end
  for (int j = 0; j < M; j++) {
    Psi[0][j] = F[0][j] + (RE/(3*DX)) * (F[1][j]-F[N-1][j]);
    Psi[N-1][j] = F[N-1][j] + (RE/(3*DX)) * (F[0][j]-F[N-2][j]);
  } // j end


  /* control matrix */
  STATIC_K = malloc_f2d(M, N);
  dlqr(J, Psi, MU, 1-MU, N, M, STATIC_K);


  free_2d((void **)J);
  free_2d((void **)F);
  free_2d((void **)Psi);
}


/* ========================================================================== */
/*   FUNCTION DEFINITIONS                                                     */
/* ========================================================================== */
/* [REQUIRED] internal setup */
void static_set() {
  /* pick from the available ROMs */
  switch (RT) {
    case BENNEY:
      benney_compute_K();
      break;
    case WR:
      ABORT("static WR not yet implemented");
      break;
    default :
      ABORT("invalid ROM type %d", RT);
  }
}

/* [REQUIRED] internal free */
void static_free() {
  free(STATIC_K);
}

/* [REQUIRED] steps the system forward in time given the interfacial height */
void static_step(double *h) {
  /* f = K * (h-1) */
  for (int i = 0; i < M; i++) {
    Amag[i] = 0.0;
    for (int j = 0; j < N; j++) {
      Amag[i] -= STATIC_K[i][j] * (interp(ITOX(j), h) - 1.0);
    } // j end
  } // i end
}

/* [REQUIRED] returns the estimator as a function of x */
double static_estimator(double x) {

  return 0.0;
}

#endif
