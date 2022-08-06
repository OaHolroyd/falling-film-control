#ifndef CONTROL_STATIC_H
#define CONTROL_STATIC_H

#include <math.h>
#include <lapacke.h>

#include "c-utils.h"
#include "linalg.h"
#include "control-core.h"


static double **STATIC_KPHI; /* control operator */


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
  double **K = malloc_f2d(M, P);
  for (i = 0; i < M; i++) {
    for (j = 0; j < P; j++) {
      K[i][j] = 0.0;
      for (k = 0; k < N; k++) {
        K[i][j] += Phi[k][j] * K_lqr[i][k] * (((double)N)/((double)P));
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
        KPHIk[i][j] += K[i][k] * Phi[j][k];
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
  while (1) {
    /* K * Phi */
    for (i = 0; i < M; i++) {
      for (j = 0; j < N; j++) {
        KPHIk[i][j] = 0.0;
        for (k = 0; k < P; k++) {
          KPHIk[i][j] += K[i][k] * Phi[j][k];
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
          Pk[i][j] += MU * KPHIk[k][i] * KPHIk[k][j];
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
    iv = 1/MU; // scale factor
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

    // solve K1 = W3\K1
    int *ipiv = malloc(P*sizeof(int));
    LAPACKE_dgesv(LAPACK_ROW_MAJOR, P, M, *W3, P, ipiv, *K1, P);


    /* update K */
    const double a = 0.005; // TODO: set optimally
    for (i = 0; i < M; i++) { // pre-solve K = W1 * W2
      for (j = 0; j < P; j++) {
        K[i][j] = (1-a)*K[i][j] + a*K1[i][j];
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
        STATIC_KPHI[i][j] += K[i][k] * Phi[j][k];
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
}


/* ========================================================================== */
/*   FUNCTION DEFINITIONS                                                     */
/* ========================================================================== */
/* [REQUIRED] internal setup */
void static_set(void) {
  STATIC_KPHI = malloc_f2d(M, N);

  /* pick from the available ROMs */
  switch (RT) {
    case BENNEY:
      static_benney_compute_KPHI(STATIC_KPHI);
      break;
    case WR:
      ABORT("WR not yet implemented");
      break;
    default :
      ABORT("invalid ROM type %d", RT);
  }
}

/* [REQUIRED] internal free */
void static_free(void) {

  free(STATIC_KPHI);
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


#endif
