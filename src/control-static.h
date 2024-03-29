#ifndef CONTROL_STATIC_H
#define CONTROL_STATIC_H

#include <math.h>
#include <lapacke.h>

#include "c-utils.h"
#include "linalg.h"
#include "control-core.h"

// #include "debug.h"

static double **STATIC_KPHI; /* control operator */
static double **STATIC_K;

/* ========================================================================== */
/*   AUXILIARY FUNCTION DEFINITIONS                                           */
/* ========================================================================== */
/* compute the control*observer matrix in the Benney case */
void static_benney_compute_KPHI(double **static_kphi) {
  int i, j, k;

  /* compute J, Psi, Phi, and K0 */
  /* Jacobian */
  double **J = malloc_f2d(N, N);
  benney_jacobian(J);

  /* actuator matrix */
  double **Psi = malloc_f2d(N, M);
  benney_actuator(Psi);

  /* observer matrix (actually the transpose) */
  double **Phi = malloc_f2d(N, P);
  benney_observer(Phi);

  /* full control matrix (note alternative cost preferences) */
  double **K_lqr = malloc_f2d(M, N);
  dlqr(J, Psi, MU, 1-MU, N, M, K_lqr);


  /* initial guess for restricted matrix */
  for (i = 0; i < M; i++) {
    for (j = 0; j < P; j++) {
      STATIC_K[i][j] = 0.0;
      for (k = 0; k < N; k++) {
        STATIC_K[i][j] += Phi[k][j] * K_lqr[i][k] * (((double)N)/((double)P));
      } // k end
    } // j end
  } // i end

  /* set up matrices for iteration */
  double **Ak = malloc_f2d(N, N);
  double **KPHIk = malloc_f2d(M, N);


  /* resulting K*Phi */
  for (i = 0; i < M; i++) {
    for (j = 0; j < N; j++) {
      KPHIk[i][j] = 0.0;
      for (k = 0; k < P; k++) {
        KPHIk[i][j] += STATIC_K[i][k] * Phi[j][k];
      } // k end
    } // j end
  } // i end

  /* closed-loop system matrix */
  for (i = 0; i < N; i++) {
    for (j = 0; j < N; j++) {
      Ak[i][j] = J[i][j];
      for (k = 0; k < M; k++) {
        Ak[i][j] -= Psi[i][k] * KPHIk[k][j];
      } // k end
    } // j end
  } // i end


  /* check that K0 is stable */
  double lmax;
  dsr(Ak, N, &lmax);
  if (lmax >= 0.0) {
    ABORT("K0 is not stable (l = %lf)", lmax);
  }

  /* iterate to find K */
  double c = 0.0; // cost
  double dc; // cost difference
  double iv;
  double **Pk = malloc_f2d(N, N);
  double **Sk = malloc_f2d(N, N);

  double **K1 = malloc_f2d(M, P); // search direction
  double **W1 = malloc_f2d(M, N); // TODO: remove if possible
  double **W2 = malloc_f2d(N, P); // TODO: remove if possible
  double **W3 = malloc_f2d(P, P); // TODO: remove if possible
  int *ipiv = malloc(P*sizeof(int));
  while (1) {
    /* K * Phi */
    for (i = 0; i < M; i++) {
      for (j = 0; j < N; j++) {
        KPHIk[i][j] = 0.0;
        for (k = 0; k < P; k++) {
          KPHIk[i][j] += STATIC_K[i][k] * Phi[j][k];
        } // k end
      } // j end
    } // i end


    /* Ak = (J - Psi * K * Phi)' */
    for (i = 0; i < N; i++) {
      for (j = 0; j < N; j++) {
        Ak[j][i] = J[i][j];
        for (k = 0; k < M; k++) {
          Ak[j][i] -= Psi[i][k] * KPHIk[k][j];
        } // k end
      } // j end
    } // i end

    /* Q = v * (K * Phi)' * K * Phi + u * I */
    for (i = 0; i < N; i++) { //
      for (j = 0; j < N; j++) {
        Pk[i][j] = 0.0;
        for (k = 0; k < M; k++) {
          Pk[i][j] += (1-MU) * KPHIk[k][i] * KPHIk[k][j];
        } // k end
      } // j end
      Pk[i][i] += DX * MU; // TODO: check the factor of DX
    } // i end

    /* solve for Pk */
    dlyap(Ak, Pk, N);


    /* Ak = J - Psi * K * Phi */
    for (i = 0; i < N; i++) {
      for (j = 0; j < N; j++) {
        Ak[i][j] = J[i][j];
        for (k = 0; k < M; k++) {
          Ak[i][j] -= Psi[i][k] * KPHIk[k][j];
        } // k end
      } // j end
    } // i end

    /* Q = I */
    for (i = 0; i < N; i++) {
      for (j = 0; j < N; j++) {
        Sk[i][j] = 0.0;
      } // j end
      Sk[i][i] = 1.0;
    } // i end

    /* solve for Sk */
    // TODO: this can be sped up singe Sk = I
    dlyap(Ak, Sk, N);


    /* compute delta cost */
    dc = 0.5 * dtr(Pk, N) - c;
    c = c + dc;


    /* compute gain update direction */
    iv = 1/(1-MU); // scale factor
    for (i = 0; i < M; i++) { // interim W1
      for (j = 0; j < N; j++) {
        W1[i][j] = 0.0;
        for (k = 0; k < N; k++) {
          W1[i][j] += iv * Psi[k][i] * Pk[k][j];
        } // k end
      } // j end
    } // i end
    for (i = 0; i < N; i++) { // interim W2
      for (j = 0; j < P; j++) {
        W2[i][j] = 0.0;
        for (k = 0; k < N; k++) {
          W2[i][j] += Sk[i][k] * Phi[k][j];
        } // k end
      } // j end
    } // i end
    for (i = 0; i < M; i++) { // pre-solve K = W1 * W2
      for (j = 0; j < P; j++) {
        K1[i][j] = 0.0;
        for (k = 0; k < N; k++) {
          K1[i][j] += W1[i][k] * W2[k][j];
        } // k end
      } // j end
    } // i end
    for (i = 0; i < P; i++) { // interim W3 (to invert)
      for (j = 0; j < P; j++) {
        W3[i][j] = 0.0;
        for (k = 0; k < N; k++) {
          W3[i][j] += Phi[k][i] * W2[k][j];
        } // k end
      } // j end
    } // i end

    // solve W3 * K1 = K1 (ie the transpose of the problem)
    LAPACKE_dgesv(LAPACK_COL_MAJOR, P, M, *W3, P, ipiv, *K1, P);


    /* update K */
    const double a = 0.005; // TODO: set optimally
    for (i = 0; i < M; i++) { // pre-solve K = W1 * W2
      for (j = 0; j < P; j++) {
        STATIC_K[i][j] = (1-a)*STATIC_K[i][j] + a*K1[i][j]; // note the transpose
      } // j end
    } // i end


    /* finish if dC is small enough */ // TODO: how small is enough
    if (fabs(dc) < 1e-4) {
      break;
    }
  }

  /* K * Phi */
  for (i = 0; i < M; i++) {
    for (j = 0; j < N; j++) {
      STATIC_KPHI[i][j] = 0.0;
      for (k = 0; k < P; k++) {
        STATIC_KPHI[i][j] += STATIC_K[i][k] * Phi[j][k];
      } // k end
    } // j end
  } // i end


  free_2d(J);
  free_2d(Psi);
  free_2d(Phi);
  free_2d(K_lqr);
  free_2d(Ak);
  free_2d(Pk);
  free_2d(Sk);
  free_2d(K1);
  free_2d(W1);
  free_2d(W2);
  free_2d(W3);
  free(ipiv);
}

/* compute the control*observer matrix in the WR case */
// void static_wr_compute_KPHI(double **static_kphi) {
//   int i, j, k;

//   /* compute J, Psi, Phi, and K0 */
//   /* Jacobian */
//   double **J = malloc_f2d(2*N, 2*N);
//   wr_jacobian(J);

//   /* actuator matrix */
//   double **Psi = malloc_f2d(2*N, M);
//   wr_actuator(Psi);

//   /* observer matrix (actually the transpose) */
//   double **Phi = malloc_f2d(2*N, 2*P);
//   wr_observer(Phi);

//   /* full control matrix (note alternative cost preferences) */
//   double **K_lqr = malloc_f2d(M, 2*N);
//   dlqr(J, Psi, MU, 1-MU, 2*N, M, K_lqr);


//   /* initial guess for restricted matrix */
//   for (i = 0; i < M; i++) {
//     for (j = 0; j < P; j++) {
//       STATIC_K[i][j] = 0.0;
//       STATIC_K[i][P+j] = 0.0;
//       for (k = 0; k < N; k++) {
//         STATIC_K[i][j] += Phi[k][j] * K_lqr[i][k] * (((double)N)/((double)P));
//         STATIC_K[i][P+j] += Phi[N+k][P+j] * K_lqr[i][N+k] * (((double)N)/((double)P));
//       } // k end
//     } // j end
//   } // i end

//   /* set up matrices for iteration */
//   double **Ak = malloc_f2d(2*N, 2*N);
//   double **KPHIk = malloc_f2d(M, 2*N);


//   /* resulting K*Phi */
//   for (i = 0; i < M; i++) {
//     for (j = 0; j < 2*N; j++) {
//       KPHIk[i][j] = 0.0;
//       for (k = 0; k < 2*P; k++) {
//         KPHIk[i][j] += STATIC_K[i][k] * Phi[j][k];
//       } // k end
//     } // j end
//   } // i end

//   /* closed-loop system matrix */
//   for (i = 0; i < 2*N; i++) {
//     for (j = 0; j < 2*N; j++) {
//       Ak[i][j] = J[i][j];
//       for (k = 0; k < M; k++) {
//         Ak[i][j] -= Psi[i][k] * KPHIk[k][j];
//       } // k end
//     } // j end
//   } // i end


//   /* check that K0 is stable */
//   double lmax;
//   dsr(Ak, N, &lmax);
//   if (lmax >= 0.0) {
//     ABORT("K0 is not stable (l = %lf)", lmax);
//   } else {
//     fprintf(stderr, "K0 is stable (l = %lf)", lmax);
//   }

//   /* iterate to find K */
//   double c = 0.0; // cost
//   double dc; // cost difference
//   double iv;
//   double **Pk = malloc_f2d(2*N, 2*N);
//   double **Sk = malloc_f2d(2*N, 2*N);

//   double **K1 = malloc_f2d(M, 2*P); // search direction
//   double **W1 = malloc_f2d(M, 2*N); // TODO: remove if possible
//   double **W2 = malloc_f2d(2*N, 2*P); // TODO: remove if possible
//   double **W3 = malloc_f2d(2*P, 2*P); // TODO: remove if possible
//   int *ipiv = malloc(2*P*sizeof(int));
//   while (1) {
//     /* K * Phi */
//     for (i = 0; i < M; i++) {
//       for (j = 0; j < 2*N; j++) {
//         KPHIk[i][j] = 0.0;
//         for (k = 0; k < 2*P; k++) {
//           KPHIk[i][j] += STATIC_K[i][k] * Phi[j][k];
//         } // k end
//       } // j end
//     } // i end


//     /* Ak = (J - Psi * K * Phi)' */
//     for (i = 0; i < 2*N; i++) {
//       for (j = 0; j < 2*N; j++) {
//         Ak[j][i] = J[i][j];
//         for (k = 0; k < M; k++) {
//           Ak[j][i] -= Psi[i][k] * KPHIk[k][j];
//         } // k end
//       } // j end
//     } // i end

//     /* Q = v * (K * Phi)' * K * Phi + u * I */
//     for (i = 0; i < 2*N; i++) { //
//       for (j = 0; j < 2*N; j++) {
//         Pk[i][j] = 0.0;
//         for (k = 0; k < M; k++) {
//           Pk[i][j] += (1-MU) * KPHIk[k][i] * KPHIk[k][j];
//         } // k end
//       } // j end
//       Pk[i][i] += DX * MU; // TODO: check the factor of DX
//     } // i end

//     /* solve for Pk */
//     dlyap(Ak, Pk, 2*N);


//     /* Ak = J - Psi * K * Phi */
//     for (i = 0; i < 2*N; i++) {
//       for (j = 0; j < 2*N; j++) {
//         Ak[i][j] = J[i][j];
//         for (k = 0; k < M; k++) {
//           Ak[i][j] -= Psi[i][k] * KPHIk[k][j];
//         } // k end
//       } // j end
//     } // i end

//     /* Q = I */
//     for (i = 0; i < 2*N; i++) {
//       for (j = 0; j < 2*N; j++) {
//         Sk[i][j] = 0.0;
//       } // j end
//       Sk[i][i] = 1.0;
//     } // i end

//     /* solve for Sk */
//     // TODO: this can be sped up singe Sk = I
//     dlyap(Ak, Sk, 2*N);


//     /* compute delta cost */
//     dc = 0.5 * dtr(Pk, 2*N) - c;
//     c = c + dc;
//     fprintf(stderr, "dc: %lf\n", dc);


//     /* compute gain update direction */
//     iv = 1/(1-MU); // scale factor
//     for (i = 0; i < M; i++) { // interim W1
//       for (j = 0; j < 2*N; j++) {
//         W1[i][j] = 0.0;
//         for (k = 0; k < 2*N; k++) {
//           W1[i][j] += iv * Psi[k][i] * Pk[k][j];
//         } // k end
//       } // j end
//     } // i end
//     for (i = 0; i < 2*N; i++) { // interim W2
//       for (j = 0; j < 2*P; j++) {
//         W2[i][j] = 0.0;
//         for (k = 0; k < 2*N; k++) {
//           W2[i][j] += Sk[i][k] * Phi[k][j];
//         } // k end
//       } // j end
//     } // i end
//     for (i = 0; i < M; i++) { // pre-solve K = W1 * W2
//       for (j = 0; j < 2*P; j++) {
//         K1[i][j] = 0.0;
//         for (k = 0; k < 2*N; k++) {
//           K1[i][j] += W1[i][k] * W2[k][j];
//         } // k end
//       } // j end
//     } // i end
//     for (i = 0; i < 2*P; i++) { // interim W3 (to invert)
//       for (j = 0; j < 2*P; j++) {
//         W3[i][j] = 0.0;
//         for (k = 0; k < 2*N; k++) {
//           W3[i][j] += Phi[k][i] * W2[k][j];
//         } // k end
//       } // j end
//     } // i end

//     // solve W3 * K1 = K1 (ie the transpose of the problem)
//     LAPACKE_dgesv(LAPACK_COL_MAJOR, 2*P, M, *W3, 2*P, ipiv, *K1, 2*P);


//     /* update K */
//     const double a = 0.005; // TODO: set optimally
//     for (i = 0; i < M; i++) { // pre-solve K = W1 * W2
//       for (j = 0; j < 2*P; j++) {
//         STATIC_K[i][j] = (1-a)*STATIC_K[i][j] + a*K1[i][j]; // note the transpose
//       } // j end
//     } // i end


//     /* finish if dC is small enough */ // TODO: how small is enough
//     if (fabs(dc) < 1e-4) {
//       break;
//     }
//   }

//   /* K * Phi */
//   for (i = 0; i < M; i++) {
//     for (j = 0; j < 2*N; j++) {
//       KPHIk[i][j] = 0.0;
//       for (k = 0; k < 2*P; k++) {
//         KPHIk[i][j] += STATIC_K[i][k] * Phi[j][k];
//       } // k end
//     } // j end
//   } // i end

//   /* apply flux approximation */
//   for (i = 0; i < M; i++) {
//     for (j = 0; j < N; j++) {
//       STATIC_KPHI[i][j] = KPHIk[i][j] + 2/3.0 * KPHIk[i][N+j];
//     } // j end
//   } // i end


//   free_2d(J);
//   free_2d(Psi);
//   free_2d(Phi);
//   free_2d(K_lqr);
//   free_2d(Ak);
//   free_2d(Pk);
//   free_2d(Sk);
//   free_2d(K1);
//   free_2d(W1);
//   free_2d(W2);
//   free_2d(W3);
//   free(ipiv);
// }


/* ========================================================================== */
/*   FUNCTION DEFINITIONS                                                     */
/* ========================================================================== */
/* [REQUIRED] internal setup */
void static_set(void) {
  STATIC_KPHI = malloc_f2d(M, N);
  STATIC_K = malloc_f2d(M, P);

  /* pick from the available ROMs */
  switch (RT) {
    case BENNEY:
      static_benney_compute_KPHI(STATIC_KPHI);
      break;
    case WR:
      // static_wr_compute_KPHI(STATIC_KPHI);
      break;
    default :
      ABORT("invalid ROM type %d", RT);
  }
}

/* [REQUIRED] internal free */
void static_free(void) {
  free(STATIC_KPHI);
  free(STATIC_K);
}

/* [REQUIRED] steps the system forward in time given the interfacial height */
void static_step(double dt, double *h) {
  /* f = K * Phi * (h-1) */
  for (int i = 0; i < M; i++) {
    Amag[i] = 0.0;
    for (int j = 0; j < N; j++) {
      Amag[i] += STATIC_KPHI[i][j] * (interp(ITOX(j), h) - 1.0);
    } // j end
  } // i end
}

/* [REQUIRED] returns the estimator as a function of x */
double static_estimator(double x) {

  return 0.0;
}

/* [REQUIRED] outputs the internal matrices */
void static_output(void) {

  output_d2d("out/K.dat", STATIC_K, M, P);
}

/* [REQUIRED] generates the control matrix CM = F*K*Phi */
void static_matrix(double **CM) {
  /* forcing matrix */
  double **F = malloc_f2d(N, M);
  forcing_matrix(F);

  for (int i = 0; i < N; i++) {
    for (int j = 0; j < N; j++) {
      CM[i][j] = 0.0;
      for (int k = 0; k < M; k++) {
        CM[i][j] += F[i][k]*STATIC_KPHI[k][j];
      } // k end
    } // j end
  } // i end

  free_2d(F);
}

#endif
