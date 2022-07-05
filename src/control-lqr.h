#ifndef CONTROL_LQR_H
#define CONTROL_LQR_H

#include <math.h>

#include "parallel.h"
#include "film-utils.h"
#include "params.h"
#include "control.h"

#include <lapacke.h> // TODO: allow for MKL


double **C_K; /* control operator */


/* ========================================================================== */
/*   AUXILIARY FUNCTION DEFINITIONS                                           */
/* ========================================================================== */
void fill_H(double **H) {
  /* forcing matrix (TODO: see if rotating this makes it faster) */
  double **F = malloc_f2d(N, C_M);
  for (int i = 0; i < N; i++) {
    for (int j = 0; j < C_M; j++) {
      F[i][j] = actuator(DX*(i+0.5)-C_loc[j]);
    } // j end
  } // i end

  /* actuator matrix (which is equal to the forcing matrix in the simplified
     Benney case) */
  double **psi = malloc_f2d(N, C_M);
  for (int i = 1; i < N-1; i++) {
    for (int j = 0; j < C_M; j++) {
      psi[i][j] = F[i][j] + (RE/(3*DX)) * (F[i+1][j]-F[i-1][j]);
    } // j end
  } // i end
  for (int j = 0; j < C_M; j++) {
    psi[0][j] = F[0][j] + (RE/(3*DX)) * (F[1][j]-F[N-1][j]);
    psi[N-1][j] = F[N-1][j] + (RE/(3*DX)) * (F[0][j]-F[N-2][j]);
  } // j end


  /* top left (J) */
  double c0 = -1/(3*CA) * (1/(DX*DX*DX*DX));
  double c1 = 1/DX + (2/tan(THETA)/3 - 8*RE/15) * (1/(DX*DX)) + 1/(3*CA) * (4/(DX*DX*DX*DX));
  double c2 = (2/tan(THETA)/3 - 8*RE/15) * (-2/(DX*DX)) - 1/(3*CA) * (6/(DX*DX*DX*DX));
  double c3 = -1/DX + (2/tan(THETA)/3 - 8*RE/15) * (1/(DX*DX)) + 1/(3*CA) * (4/(DX*DX*DX*DX));
  double c4 = -1/(3*CA) * (1/(DX*DX*DX*DX));
  for (int i = 2; i < N-2; i++) {
    H[i][i-2] = c0;
    H[i][i-1] = c1;
    H[i][i] = c2;
    H[i][i+1] = c3;
    H[i][i+2] = c4;
  } // i end

  H[0][N-2] = c0;
  H[0][N-1] = c1;
  H[0][0] = c2;
  H[0][1] = c3;
  H[0][2] = c4;

  H[1][N-1] = c0;
  H[1][0] = c1;
  H[1][1] = c2;
  H[1][1+1] = c3;
  H[1][1+2] = c4;

  H[N-2][N-4] = c0;
  H[N-2][N-3] = c1;
  H[N-2][N-2] = c2;
  H[N-2][N-1] = c3;
  H[N-2][0] = c4;

  H[N-1][N-3] = c0;
  H[N-1][N-2] = c1;
  H[N-1][N-1] = c2;
  H[N-1][0] = c3;
  H[N-1][1] = c4;


  /* bottom right (-J^T) */
  c4 = 1/(3*CA) * (1/(DX*DX*DX*DX));
  c3 = -1/DX - (2/tan(THETA)/3 - 8*RE/15) * (1/(DX*DX)) - 1/(3*CA) * (4/(DX*DX*DX*DX));
  c2 = -(2/tan(THETA)/3 - 8*RE/15) * (-2/(DX*DX)) + 1/(3*CA) * (6/(DX*DX*DX*DX));
  c1 = 1/DX - (2/tan(THETA)/3 - 8*RE/15) * (1/(DX*DX)) - 1/(3*CA) * (4/(DX*DX*DX*DX));
  c0 = 1/(3*CA) * (1/(DX*DX*DX*DX));
  for (int i = 2; i < N-2; i++) {
    H[N+i][N+i-2] = c0;
    H[N+i][N+i-1] = c1;
    H[N+i][N+i] = c2;
    H[N+i][N+i+1] = c3;
    H[N+i][N+i+2] = c4;
  } // i end

  H[N+0][N+N-2] = c0;
  H[N+0][N+N-1] = c1;
  H[N+0][N+0] = c2;
  H[N+0][N+1] = c3;
  H[N+0][N+2] = c4;

  H[N+1][N+N-1] = c0;
  H[N+1][N+0] = c1;
  H[N+1][N+1] = c2;
  H[N+1][N+1+1] = c3;
  H[N+1][N+1+2] = c4;

  H[N+N-2][N+N-4] = c0;
  H[N+N-2][N+N-3] = c1;
  H[N+N-2][N+N-2] = c2;
  H[N+N-2][N+N-1] = c3;
  H[N+N-2][N+0] = c4;

  H[N+N-1][N+N-3] = c0;
  H[N+N-1][N+N-2] = c1;
  H[N+N-1][N+N-1] = c2;
  H[N+N-1][N+0] = c3;
  H[N+N-1][N+1] = c4;


  /* bottom left (-U) */
  for (int i = 0; i < N; i++) {
    H[N+i][i] = -DX*C_MU;
  } // i end


  /* top right (-PHI * V^-1 * PSI^T) */
  c0 = 1.0/(C_MU-1.0);
  for (int i = 0; i < N; i++) {
    for (int j = i; j < N; j++) {
      H[i][N+j] = 0.0;
      for (int k = 0; k < C_M; k++) {
        H[i][N+j] += psi[i][k] * psi[j][k];
      } // k end
      H[i][N+j] *= c0;
      H[j][N+i] = H[i][N+j]; // use symmetry of PHI * PSI^T
    } // j end
  } // i end

  free_2d((void **)psi);
}

/* ========================================================================== */
/*   FUNCTION DEFINITIONS                                                     */
/* ========================================================================== */
/* [REQUIRED] sets the control variables (call after set_params) */
void set_Cparams() {
  internal_set_Cparams();

  /* control operator (BENNEY SYSTEM) */
  C_K = malloc_f2d(C_M, N);

  /* interim arrays */
  double **H = malloc_f2d(2*N, 2*N);
  fill_H(H);

  /* compute eigenvalues/vectors */
  double *wr = malloc(2*N*sizeof(double));
  double *wi = malloc(2*N*sizeof(double));
  double *vr = malloc(2*N*2*N*sizeof(double));
  int info = LAPACKE_dgeev(LAPACK_ROW_MAJOR, 'N', 'V', 2*N, *H, 2*N, wr, wi, NULL, 1, vr, 2*N);

  FILE *fp = fopen("out/H.dat", "w");
  for (int i = 0; i < 2*N; i++) {
    for (int j = 0; j < 2*N; j++) {
      fprintf(fp, "%g ", H[i][j]);
    } // j end
    fprintf(fp, "\n");
  } // i end
  fclose(fp);


  /* free workspace */
  free_2d((void **)H);
  free(wr);
  free(wi);
  free(vr);

  fprintf(stderr, "set K\n");
  abort();
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
    C_mag[i] = interfacial_height(C_loc[i] - C_PHI) - 1;
  } // i end
}


#endif
