#ifndef CONTROL_LQR_H
#define CONTROL_LQR_H

#include <math.h>

#include "c-utils.h"
#include "parallel.h"
#include "film-utils.h"
#include "params.h"
#include "control.h"
#include "lqr.h"


COMPLEX **C_K; // control operator
COMPLEX **C_A; // spectral dynamics A
COMPLEX **C_B; // spectral dynamics B
COMPLEX *C_z; // estimator
double *C_Oloc; // observer locations

COMPLEX **C_D; // inversion matrix
COMPLEX *C_z0; // previous timestep
int *C_ipiv; // pivot array


/* ========================================================================== */
/*   AUXILIARY FUNCTION DEFINITIONS                                           */
/* ========================================================================== */
void update_estimator_benney() {
  int i, j, k;

  /* store current estimator */
  for (i = 0; i < C_M; i++) {
    C_z0[i] = C_z[i];
  } // i end


  /* set implicit matrix */
  for (i = 0; i < C_M; i++) {
    for (j = 0; j < C_M; j++) {
      C_D[i][j] = -0.5*dt*C_A[i][j];
    } // j end
  } // i end

  for (i = 0; i < C_M; i++) {
    C_D[i][i] += 1.0;
  } // i end


  /* compute the explicit half */
  for (i = 0; i < C_M; i++) {
    for (j = 0; j < C_M; j++) {
      C_z[i] += 0.5*dt*C_A[i][j]*C_z0[j];
      // C_z[i] += dt*C_A[i][j]*C_z0[j];
    } // j end
  } // i end


  /* compute the explicit half */
  for (i = 0; i < C_M; i++) {
    for (j = 0; j < N; j++) {
      C_z[i] += dt*C_B[i][j]*(interfacial_height(DX*(j+0.5)) - 1);
    } // j end
  } // i end


  /* solve the implicit half */
  LQR_zgesv(C_M, 1, *C_D, C_M, C_ipiv, C_z, C_M);
}


/* ========================================================================== */
/*   FUNCTION DEFINITIONS                                                     */
/* ========================================================================== */
/* [REQUIRED] sets the control variables (call after set_params) */
void set_Cparams() {
  internal_set_Cparams();

  int i, j, k;
  FILE *fp; // TODO: remove outputs


  /* observer locations */
  C_Oloc = malloc(C_P*sizeof(double));
  for (i = 0; i < C_P; i++) {
    C_Oloc[i] = (i+0.5)*LX/C_P;
  } // i end


  /* compute (real) observer matrix (actually the transpose) */
  double **r_Phi = malloc_f2d(N, C_P);
  for (i = 0; i < N; i++) {
    for (j = 0; j < C_P; j++) {
      r_Phi[i][j] = actuator(DX*(i+0.5)-C_Oloc[j]);
    } // j end
  } // i end

  // fp = fopen("out/r_Phi.dat", "w");
  // fprintf(stderr, "done r_Phi\n");
  // for (i = 0; i < N; i++) {
  //   for (j = 0; j < C_P; j++) {
  //     fprintf(fp, "%lf ", r_Phi[i][j]);
  //   } // j end
  //   fprintf(fp, "\n");
  // } // i end
  // fclose(fp);


  /* compute controllable wavenumbers */
  double *K = malloc(C_M*sizeof(double));
  double dk = 2 * M_PI / LX;
  j = 1;
  for (i = 0; i < C_M; i++) {
    j *= -1;
    k = j*(i+1)/2; // uses integer rounding
    K[i] = dk*k;
  } // i end


  /* compute spectral Jacobian */
  COMPLEX **J = malloc_z2d(C_M, C_M);
  for (i = 0; i < C_M; i++) {
    for (j = 0; j < C_M; j++) {
      J[i][j] = 0.0;
    } // j end
    J[i][i] = -2*I*K[i] - (2.0/3.0/tan(THETA) - 8.0*RE/15.0)*K[i]*K[i] - 1.0/3.0/CA*K[i]*K[i]*K[i]*K[i];
  } // i end

  // fp = fopen("out/Jr.dat", "w");
  // fprintf(stderr, "done Jr\n");
  // for (i = 0; i < C_M; i++) {
  //   for (j = 0; j < C_M; j++) {
  //     fprintf(fp, "%lf ", creal(J[i][j]));
  //   } // j end
  //   fprintf(fp, "\n");
  // } // i end
  // fclose(fp);

  // fp = fopen("out/Ji.dat", "w");
  // for (i = 0; i < C_M; i++) {
  //   for (j = 0; j < C_M; j++) {
  //     fprintf(fp, "%lf ", cimag(J[i][j]));
  //   } // j end
  //   fprintf(fp, "\n");
  // } // i end
  // fclose(fp);


  /* compute spectral Psi */
  COMPLEX **Psi = malloc_z2d(C_M, C_M);
  double xk;
  for (i = 0; i < C_M; i++) {
    for (j = 0; j < C_M; j++) {
      Psi[i][j] = 0.0;

      /* manually compute the DFT via dot product */
      for (k = 0; k < N; k++) {
        xk = DX*(k+0.5);
        Psi[i][j] += actuator(xk-C_loc[j]) * (cos(K[i]*xk) - I*sin(K[i]*xk));
      } // k end

      /* (1 + 2RE/3 dx) component */
      Psi[i][j] *= 1.0 + 2.0*RE/3.0*I*K[i];
    } // j end
  } // i end

  // fp = fopen("out/Psir.dat", "w");
  // fprintf(stderr, "done Psir\n");
  // for (i = 0; i < C_M; i++) {
  //   for (j = 0; j < C_M; j++) {
  //     fprintf(fp, "%lf ", creal(Psi[i][j]));
  //   } // j end
  //   fprintf(fp, "\n");
  // } // i end
  // fclose(fp);

  // fp = fopen("out/Psii.dat", "w");
  // for (i = 0; i < C_M; i++) {
  //   for (j = 0; j < C_M; j++) {
  //     fprintf(fp, "%lf ", cimag(Psi[i][j]));
  //   } // j end
  //   fprintf(fp, "\n");
  // } // i end
  // fclose(fp);


  /* compute spectral Phi (this is actually Phi') */
  COMPLEX **Phi = malloc_z2d(C_M, C_P);
  for (i = 0; i < C_M; i++) {
    for (j = 0; j < C_P; j++) {
      Phi[i][j] = 0.0;

      /* manually compute the DFT via dot product */
      for (k = 0; k < N; k++) {
        xk = DX*(k+0.5);
        Phi[i][j] += actuator(xk-C_Oloc[j]) * (cos(K[i]*xk) - I*sin(K[i]*xk));
      } // k end
    } // j end
  } // i end

  // fp = fopen("out/Phir.dat", "w");
  // fprintf(stderr, "done Phir\n");
  // for (i = 0; i < C_M; i++) {
  //   for (j = 0; j < C_P; j++) {
  //     fprintf(fp, "%lf ", creal(Phi[i][j]));
  //   } // j end
  //   fprintf(fp, "\n");
  // } // i end
  // fclose(fp);

  // fp = fopen("out/Phii.dat", "w");
  // for (i = 0; i < C_M; i++) {
  //   for (j = 0; j < C_P; j++) {
  //     fprintf(fp, "%lf ", cimag(Phi[i][j]));
  //   } // j end
  //   fprintf(fp, "\n");
  // } // i end
  // fclose(fp);


  /* compute spectral K */
  C_K = malloc_z2d(C_M, C_M);
  zlqr(J, Psi, C_MU*LX/(C_M*C_M), (1-C_MU), C_M, C_M, C_K);

  // fp = fopen("out/Kr.dat", "w");
  // fprintf(stderr, "done Kr\n");
  // for (i = 0; i < C_M; i++) {
  //   for (j = 0; j < C_M; j++) {
  //     fprintf(fp, "%lf ", creal(C_K[i][j]));
  //   } // j end
  //   fprintf(fp, "\n");
  // } // i end
  // fclose(fp);

  // fp = fopen("out/Ki.dat", "w");
  // for (i = 0; i < C_M; i++) {
  //   for (j = 0; j < C_M; j++) {
  //     fprintf(fp, "%lf ", cimag(C_K[i][j]));
  //   } // j end
  //   fprintf(fp, "\n");
  // } // i end
  // fclose(fp);


  /* compute L (this is actually L.') */
  COMPLEX **L = malloc_z2d(C_P, C_M);
  zlqr(J, Phi, C_MU*LX/(C_M*C_M), (1-C_MU), C_M, C_P, L);

  // fp = fopen("out/Lr.dat", "w");
  // fprintf(stderr, "done Lr\n");
  // for (i = 0; i < C_P; i++) {
  //   for (j = 0; j < C_M; j++) {
  //     fprintf(fp, "%lf ", creal(L[i][j]));
  //   } // j end
  //   fprintf(fp, "\n");
  // } // i end
  // fclose(fp);

  // fp = fopen("out/Li.dat", "w");
  // for (i = 0; i < C_P; i++) {
  //   for (j = 0; j < C_M; j++) {
  //     fprintf(fp, "%lf ", cimag(L[i][j]));
  //   } // j end
  //   fprintf(fp, "\n");
  // } // i end
  // fclose(fp);


  /* compute B */
  C_B = malloc_z2d(C_M, N);
  for (i = 0; i < C_M; i++) {
    for (j = 0; j < N; j++) {
      C_B[i][j] = 0.0;
      for (k = 0; k < C_P; k++) {
        C_B[i][j] -= L[k][i] * r_Phi[j][k];
      } // k end
    } // j end
  } // i end

  // fp = fopen("out/Br.dat", "w");
  // fprintf(stderr, "done Br\n");
  // for (i = 0; i < C_M; i++) {
  //   for (j = 0; j < N; j++) {
  //     fprintf(fp, "%lf ", creal(C_B[i][j]));
  //   } // j end
  //   fprintf(fp, "\n");
  // } // i end
  // fclose(fp);

  // fp = fopen("out/Bi.dat", "w");
  // for (i = 0; i < C_M; i++) {
  //   for (j = 0; j < N; j++) {
  //     fprintf(fp, "%lf ", cimag(C_B[i][j]));
  //   } // j end
  //   fprintf(fp, "\n");
  // } // i end
  // fclose(fp);


  /* compute A */
  C_A = malloc_z2d(C_M, C_M);
  for (i = 0; i < C_M; i++) {
    for (j = 0; j < C_M; j++) {
      C_A[i][j] = J[i][j];
      for (k = 0; k < C_M; k++) {
        C_A[i][j] += Psi[i][k] * C_K[k][j];
      } // k end
      for (k = 0; k < C_P; k++) {
        C_A[i][j] += L[k][i] * Phi[j][k];
      } // k end
    } // j end
  } // i end

  // fp = fopen("out/Ar.dat", "w");
  // fprintf(stderr, "done Ar\n");
  // for (i = 0; i < C_M; i++) {
  //   for (j = 0; j < C_M; j++) {
  //     fprintf(fp, "%lf ", creal(C_A[i][j]));
  //   } // j end
  //   fprintf(fp, "\n");
  // } // i end
  // fclose(fp);

  // fp = fopen("out/Ai.dat", "w");
  // for (i = 0; i < C_M; i++) {
  //   for (j = 0; j < C_M; j++) {
  //     fprintf(fp, "%lf ", cimag(C_A[i][j]));
  //   } // j end
  //   fprintf(fp, "\n");
  // } // i end
  // fclose(fp);


  /* allocate D, z, and z0 */
  C_D = malloc_z2d(C_M, C_M);
  C_z = malloc(C_M * sizeof(COMPLEX));
  C_z0 = malloc(C_M * sizeof(COMPLEX));
  C_ipiv = malloc(C_M*sizeof(int));

  for (i = 0; i < C_M; i++) {
    C_z[i] = 0.0;
  } // i end


  /* free J, Psi, Phi, L */
  free_2d(J);
  free_2d(Psi);
  free_2d(Phi);
  free_2d(r_Phi);
  free_2d(L);
  free(K);

  // ABORT("normal end");
}

/* [REQUIRED] frees control variables */
void control_free() {
  internal_control_free();

  free_2d(C_K);
  free_2d(C_A);
  free_2d(C_B);
  free_2d(C_D);
  free(C_z);
  free(C_z0);
  free(C_Oloc);
}

/* [REQUIRED] sets the control magnitudes and returns the cost */
void control_set_magnitudes() {
  update_estimator_benney();

  /* basic paired observer-actuators */
  for (int i = 0; i < C_M; i++) {
    C_mag[i] = 0.0;
    for (int j = 0; j < C_M; j++) {
      C_mag[i] -= creal(C_K[i][j] * C_z[j]);
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
  COMPLEX e = 0.0;
  double dk = 2 * M_PI / LX;
  int j, k;
  j = 1;
  for (int i = 0; i < C_M; i++) {
    j *= -1;
    k = j*(i+1)/2; // uses integer rounding
    e += C_z[i]*(cos(dk*k*x) - I*sin(dk*k*x));
  } // i end

  return creal(e);
}

#endif
