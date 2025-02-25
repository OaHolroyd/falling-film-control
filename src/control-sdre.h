#ifndef CONTROL_SDRE_H
#define CONTROL_SDRE_H

#include <math.h>
#include <string.h>

#include "c-utils.h"
#include "linalg.h"
#include "control-core.h"


static double **SDRE_K1; // first-order gain operator
static double ***SDRE_K2; // basis for second-order gain correction


/* ========================================================================== */
/*   AUXILIARY FUNCTION DEFINITIONS                                           */
/* ========================================================================== */
/**
 * Wrap an index using periodic BCs
 * @param  k index to wrap
 * @return   wrapped index
 */
int wrap(int k) {
  if (k >= N) {
    return wrap(k-N);
  } else if (k < 0) {
    return wrap(k+N);
  } else {
    return k;
  }
}

/* construct the constant part of the higher order corrective term:
    A(h) = A0 + A1(h),

             N
    A1(h) = SUM fk(h) * dAk,
            k=1

    fk(h) = hk.
 */
void construct_dAi_benney(double **dAk, int k) {
  // most of the matrix will b zero
  memset(*dAk, 0, N*N*sizeof(double));

  const double BETA = 1.0 / tan(THETA);

  dAk[wrap(k-1)][wrap(k-3)] = - 1.0/CA * -0.25/DX/DX/DX/DX;
  dAk[wrap(k-1)][wrap(k-2)] = (-16.0/5.0 * RE + 2.0 * BETA) * -0.25/DX/DX - 1.0/CA * 0.5/DX/DX/DX/DX;
  dAk[wrap(k-1)][wrap(k+0)] = (-16.0/5.0 * RE + 2.0 * BETA) * 0.25/DX/DX - 1.0/CA * -0.5/DX/DX/DX/DX;
  dAk[wrap(k-1)][wrap(k+1)] = - 1.0/CA * 0.25/DX/DX/DX/DX;

  dAk[wrap(k)][wrap(k-2)] = - 1/CA * 1/DX/DX/DX/DX;
  dAk[wrap(k)][wrap(k-1)] = (-16.0/5.0 * RE + 2.0 * BETA) * 1.0/DX/DX - 4.0 * -0.5/DX - 1.0/CA * -4.0/DX/DX/DX/DX;
  dAk[wrap(k)][wrap(k+0)] = (-16.0/5.0 * RE + 2.0 * BETA) * -2.0/DX/DX - 1.0/CA * 6.0/DX/DX/DX/DX;
  dAk[wrap(k)][wrap(k+1)] = (-16.0/5.0 * RE + 2.0 * BETA) * 1.0/DX/DX - 4.0 * 0.5/DX - 1.0/CA * -4.0/DX/DX/DX/DX;
  dAk[wrap(k)][wrap(k+2)] =  - 1.0/CA * 1.0/DX/DX/DX/DX;

  dAk[wrap(k+1)][wrap(k-1)] = - 1.0/CA * 0.25/DX/DX/DX/DX;
  dAk[wrap(k+1)][wrap(k+0)] = (-16.0/5.0 * RE + 2.0 * BETA) * 0.25/DX/DX - 1.0/CA * -0.5/DX/DX/DX/DX;
  dAk[wrap(k+1)][wrap(k+2)] = (-16.0/5.0 * RE + 2.0 * BETA) * -0.25/DX/DX - 1.0/CA * 0.5/DX/DX/DX/DX;
  dAk[wrap(k+1)][wrap(k+3)] = - 1.0/CA * -0.25/DX/DX/DX/DX;
}


void construct_dAi_wr(double **dAk, int k) {
  // most of the matrix will b zero
  memset(*dAk, 0, 4*N*N*sizeof(double));

  if (k < N) {
    // bottom left quadrant (h->q)
    dAk[N+k][wrap(k-2)] = 5.0/4.0/RE/CA * -0.5/DX/DX/DX;
    dAk[N+k][wrap(k-1)] = -5.0/2.0/RE * -0.5/DX + 5.0/4.0/RE/CA * 1.0/DX/DX/DX;
    dAk[N+k][wrap(k+0)] = 5.0/2.0/RE * 1.0;
    dAk[N+k][wrap(k+1)] = -5.0/2.0/RE * 0.5/DX + 5.0/4.0/RE/CA * -1.0/DX/DX/DX;
    dAk[N+k][wrap(k+2)] = 5.0/4.0/RE/CA * 0.5/DX/DX/DX;

    // bottom right quadrant (q->q)
    dAk[N+k][N+wrap(k-1)] = -17.0/21.0 * -0.5/DX;
    dAk[N+k][N+wrap(k+1)] = -17.0/21.0 * 0.5/DX;
  } else {
    k -= N;
    // bottom left quadrant (h->q)
    dAk[N+k][wrap(k-1)] = 6.0/7.0 * -0.5/DX;
    dAk[N+k][wrap(k+1)] = 6.0/7.0 * 0.5/DX;

    // bottom right quadrant (q->q)
    dAk[N+k][N+wrap(k-1)] = -17.0/14.0 * -0.5/DX;
    dAk[N+k][N+wrap(k+1)] = -17.0/14.0 * 0.5/DX;
  }
}


/* compute the fixed gain operator and the basis for the gain correction in the Benney case */
void sdre_benney_compute_matrices(void) {

  /* Jacobian */
  double **A = malloc_f2d(N, N);
  benney_jacobian(A);

  /* actuator matrix */
  double **B = malloc_f2d(N, M);
  benney_actuator(B);

  /* first-order gain matrix */
  dlqr(A, B, DX*MU, 1.0-MU, N, M, SDRE_K1);
  // output_d2d("out/K1.dat", SDRE_K1, M, N);

  // also need to know L0
  double complex **L0_complex = malloc_z2d(N, N);
  riccati(A, B, DX*MU, 1.0-MU, N, M, *L0_complex);
  // convert to real
  double **L0 = malloc_f2d(N, N);
  for (int i = 0; i < N; i++) {
    for (int j = 0; j < N; j++) {
      L0[i][j] = creal(L0_complex[i][j]);
    }
  }
  // output_d2d("out/L0.dat", L0, N, N);

  // compute A_lyap, the fixed term in the computation of the basis matrices
  double **A_lyap = malloc_f2d(N, N);
  for (int ii = 0; ii < N; ii++) {
    for (int jj = 0; jj < N; jj++) {
      A_lyap[ii][jj] = A[ii][jj];
      for (int kk = 0; kk < M; kk++) {
        A_lyap[ii][jj] -= B[ii][kk] * SDRE_K1[kk][jj];
      }
    }
  }
  // output_d2d("out/A_lyap.dat", A_lyap, N, N);

  // for each element of the vector h, compute a basis matrix
  double **dAi = malloc_f2d(N, N);
  double **Q_lyap = malloc_f2d(N, N);
  for (int i = 0; i < N; i++) {
    // compute ∆Ai for the ith nonlinear component of the Benney matrix
    construct_dAi_benney(dAi, i);

    char filename[128];
    // sprintf(filename, "out/A2-%d.dat", i);
    // output_d2d(filename, dAi, N, N);


    // compute Q_lyap = L0 * ∆Ai + ∆Ai * L0
    for (int ii = 0; ii < N; ii++) {
      for (int jj = 0; jj < N; jj++) {
        Q_lyap[ii][jj] = 0.0;
        for (int kk = 0; kk < N; kk++) {
          Q_lyap[ii][jj] += L0[ii][kk] * dAi[kk][jj];
          Q_lyap[ii][jj] += dAi[ii][kk] * L0[kk][jj];
        }
      }
    }

    for (int ii = 0; ii < N; ii++) {
      for (int jj = 0; jj < N; jj++) {
        A_lyap[ii][jj] = A[ii][jj];
        for (int kk = 0; kk < M; kk++) {
          A_lyap[ii][jj] -= B[ii][kk] * SDRE_K1[kk][jj];
        }
      }
    }

    // sprintf(filename, "out/Qlyap-%d.dat", i);
    // output_d2d(filename, Q_lyap, N, N);

    // solve the Lyapunov problem for Li
    dlyap(A_lyap, Q_lyap, N, 0);
    double **Li = Q_lyap; // just for notational convenience

    // sprintf(filename, "out/L1-%d.dat", i);
    // output_d2d(filename, Li, N, N);

    // compute gain version Ki = 1 / (1-MU) * B' * Li
    double **Ki = SDRE_K2[i];
    for (int ii = 0; ii < M; ii++) {
      for (int jj = 0; jj < N; jj++) {
        Ki[ii][jj] = 0.0;
        for (int kk = 0; kk < N; kk++) {
          Ki[ii][jj] += B[kk][ii] * Li[kk][jj];
        }
        Ki[ii][jj] *= 1.0 / (1.0 - MU);
      }
    }

    sprintf(filename, "out/K2-%d.dat", i);
    output_d2d(filename, Ki, M, N);
  }

  free_2d(A);
  free_2d(B);
  free_2d(L0_complex);
  free_2d(A_lyap);
  free_2d(Q_lyap);
  free_2d(dAi);
}


/* compute the fixed gain operator and the basis for the gain correction in the WR case */
void sdre_wr_compute_matrices(void) {
  /* Jacobian */
  double **A = malloc_f2d(2*N, 2*N);
  wr_jacobian(A);

  /* actuator matrix */
  double **B = malloc_f2d(2*N, M);
  wr_actuator(B);

  /* first-order gain matrix */
  dlqr(A, B, DX*MU, 1.0-MU, 2*N, M, SDRE_K1);
  // output_d2d("out/K1.dat", SDRE_K1, M, 2*N);

  // also need to know L0
  double complex **L0_complex = malloc_z2d(2*N, 2*N);
  riccati(A, B, DX*MU, 1.0-MU, 2*N, M, *L0_complex);
  // convert to real
  double **L0 = malloc_f2d(2*N, 2*N);
  for (int i = 0; i < 2*N; i++) {
    for (int j = 0; j < 2*N; j++) {
      L0[i][j] = creal(L0_complex[i][j]);
    }
  }
  // output_d2d("out/L0.dat", L0, 2*N, 2*N);

  // compute A_lyap, the fixed term in the computation of the basis matrices
  double **A_lyap = malloc_f2d(2*N, 2*N);
  for (int ii = 0; ii < 2*N; ii++) {
    for (int jj = 0; jj < 2*N; jj++) {
      A_lyap[ii][jj] = A[ii][jj];
      for (int kk = 0; kk < M; kk++) {
        A_lyap[ii][jj] -= B[ii][kk] * SDRE_K1[kk][jj];
      }
    }
  }
  // output_d2d("out/A_lyap.dat", A_lyap, 2*N, 2*N);

  // for each element of the vector [h, q], compute a basis matrix
  double **dAi = malloc_f2d(2*N, 2*N);
  double **Q_lyap = malloc_f2d(2*N, 2*N);
  for (int i = 0; i < 2*N; i++) {
    // compute ∆Ai for the ith nonlinear component of the WR matrix
    construct_dAi_wr(dAi, i);

    // char filename[128];
    // sprintf(filename, "out/A2-%d.dat", i);
    // output_d2d(filename, dAi, 2*N, 2*N);


    // compute Q_lyap = L0 * ∆Ai + ∆Ai * L0
    for (int ii = 0; ii < 2*N; ii++) {
      for (int jj = 0; jj < 2*N; jj++) {
        Q_lyap[ii][jj] = 0.0;
        for (int kk = 0; kk < 2*N; kk++) {
          Q_lyap[ii][jj] += L0[ii][kk] * dAi[kk][jj];
          Q_lyap[ii][jj] += dAi[ii][kk] * L0[kk][jj];
        }
      }
    }

    // recompute A_lyap
    for (int ii = 0; ii < 2*N; ii++) {
      for (int jj = 0; jj < 2*N; jj++) {
        A_lyap[ii][jj] = A[ii][jj];
        for (int kk = 0; kk < M; kk++) {
          A_lyap[ii][jj] -= B[ii][kk] * SDRE_K1[kk][jj];
        }
      }
    }

    // sprintf(filename, "out/Qlyap-%d.dat", i);
    // output_d2d(filename, Q_lyap, 2*N, 2*N);

    // solve the Lyapunov problem for Li
    dlyap(A_lyap, Q_lyap, 2*N, 0);
    double **Li = Q_lyap; // just for notational convenience

    // sprintf(filename, "out/L1-%d.dat", i);
    // output_d2d(filename, Li, 2*N, 2*N);

    // compute gain version Ki = 1 / (1-MU) * B' * Li
    double **Ki = SDRE_K2[i];
    for (int ii = 0; ii < M; ii++) {
      for (int jj = 0; jj < 2*N; jj++) {
        Ki[ii][jj] = 0.0;
        for (int kk = 0; kk < 2*N; kk++) {
          Ki[ii][jj] += B[kk][ii] * Li[kk][jj];
        }
        Ki[ii][jj] *= 1.0 / (1.0 - MU);
      }
    }

    // sprintf(filename, "out/K2-%d.dat", i);
    // output_d2d(filename, Ki, M, 2*N);
  }

  free_2d(A);
  free_2d(B);
  free_2d(L0_complex);
  free_2d(A_lyap);
  free_2d(Q_lyap);
  free_2d(dAi);
}


/* ========================================================================== */
/*   FUNCTION DEFINITIONS                                                     */
/* ========================================================================== */
/* [REQUIRED] internal setup */
void sdre_set(void) {
  /* pick from the available ROMs */
  switch (RT) {
    case BENNEY: {
        // allocate space
        SDRE_K1 = malloc_f2d(M, N);
        SDRE_K2 = malloc(N * sizeof(double**));
        for (int i = 0; i < N; i++) {
          // TODO: allocate this in a proper block rather than in a loop
          SDRE_K2[i] = malloc_f2d(M, N);
        }
        sdre_benney_compute_matrices();
        break;
      }
    case WR: {
        // allocate space
        SDRE_K1 = malloc_f2d(M, 2*N);
        SDRE_K2 = malloc(2*N * sizeof(double**));
        for (int i = 0; i < 2*N; i++) {
          // TODO: allocate this in a proper block rather than in a loop
          SDRE_K2[i] = malloc_f2d(M, 2*N);
        }

        sdre_wr_compute_matrices();
        break;
     }
    default :
      ABORT("invalid ROM type %d", RT);
  }
}

/* [REQUIRED] internal free */
void sdre_free(void) {
  free_2d(SDRE_K1);

  for (int i = 0; i < N; i++) {
    free_2d(SDRE_K2[i]);
  }

  if (RT == WR) {
    for (int i = N; i < 2*N; i++) {
      free_2d(SDRE_K2[i]);
    }
  }

  free_2d(SDRE_K2);
}

/* [REQUIRED] steps the system forward in time given the interfacial height */
void sdre_step(double dt, double *h, double *q) {
  /* f = K1 * (h-1) */
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

  // fprintf(stderr, "\n\n");
  // for (int i = 0; i < M; i++) {
  //   fprintf(stderr, "%lf ", Amag[i]);
  // }

  /* quadratic correction */
  double scale = -1.0 / sqrt(((double)(N))); // needs to be negative
  for (int k = 0; k < N; k++) {
    double fk = scale * (interp(ITOX(k), h) - 1.0);
    for (int i = 0; i < M; i++) {
      for (int j = 0; j < N; j++) {
        Amag[i] += fk * SDRE_K2[k][i][j] * (interp(ITOX(j), h) - 1.0);
      } // j end

      /* only WR uses the flux */
      if (RT == WR) {
        for (int j = 0; j < N; j++) {
          Amag[i] += fk * SDRE_K2[k][i][N+j] * (interp(ITOX(j), q) - 2.0/3.0);
        } // j end
      }
    } // i end
  }

  /* WR has more correction terms */
  if (RT == WR) {
    for (int k = 0; k < N; k++) {
      double fk = scale * (interp(ITOX(k), q) - 2.0/3.0);
      for (int i = 0; i < M; i++) {
        for (int j = 0; j < N; j++) {
          Amag[i] += fk * SDRE_K2[N+k][i][j] * (interp(ITOX(j), h) - 1.0);
        } // j end

        for (int j = 0; j < N; j++) {
          Amag[i] += fk * SDRE_K2[N+k][i][N+j] * (interp(ITOX(j), q) - 2.0/3.0);
        } // j end
      } // i end
    }
  }


  // fprintf(stderr, "\n");
  // for (int i = 0; i < M; i++) {
  //   fprintf(stderr, "%lf ", Amag[i]);
  // }
  // fprintf(stderr, "\n\n");

  // fprintf(stderr, "h = [");
  // for (int i = 0; i < N-1; i++) {
  //   fprintf(stderr, "%.16lf, ", (interp(ITOX(i), h) - 1.0));
  // }
  // fprintf(stderr, "%.16lf]';\n", (interp(ITOX(N-1), h) - 1.0));

  // exit(1);
}

/* [REQUIRED] returns the estimator as a function of x */
double sdre_estimator(double x) {

  return 0.0;
}

/* [REQUIRED] outputs the internal matrices */
void sdre_output(void) {
  // if (RT == BENNEY) {
  //   output_d2d("out/K.dat", SDRE_K1, M, N);
  // } else {
  //   output_d2d("out/K.dat", SDRE_K1, M, 2*N);
  // }
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
