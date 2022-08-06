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
void lqr_benney_compute_K(void) {
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
  dlqr(J, Psi, MU/DX, 1-MU, N, M, LQR_K);


  free_2d(J);
  free_2d(F);
  free_2d(Psi);
}

/* compute the control matrix in the weighted-residuals case */
void lqr_wr_compute_K(void) {
    /* Jacobian */
  double **J = malloc_f2d(2*N, 2*N);

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


  /* forcing matrix (TODO: see if rotating this makes it faster) */
  double **F = malloc_f2d(N, M);
  for (int i = 0; i < N; i++) {
    for (int j = 0; j < M; j++) {
      F[i][j] = actuator(ITOX(i)-Aloc[j]);
    } // j end
  } // i end


  /* actuator matrix */
  double **Psi = malloc_f2d(2*N, M);
  for (int i = 0; i < N; i++) {
    for (int j = 0; j < M; j++) {
      Psi[i][j] = F[i][j];
      Psi[i+N][j] = 2*RE/15 * F[i][j];
    } // j end
  } // i end


  /* full control matrix */
  double **K = malloc_f2d(M, 2*N);
  dlqr(J, Psi, MU, 1-MU, 2*N, M, K);

  /* apply flux approximation */
  for (int i = 0; i < M; i++) {
    for (int j = 0; j < N; j++) {
      LQR_K[i][j] = K[i][j] + 2/3 * K[i][j+N];
    } // j end
  } // i end


  free_2d(J);
  free_2d(F);
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
      lqr_benney_compute_K();
      break;
    case WR:
      lqr_wr_compute_K();
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
