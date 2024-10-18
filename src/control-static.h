#ifndef CONTROL_STATIC_H
#define CONTROL_STATIC_H

#include <stdio.h>
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
    dlyap(Ak, Pk, N, 0);


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
    dlyap(Ak, Sk, N, 0);


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
void static_wr_compute_KPHI(double **static_kphi) {
  // fprintf(stderr, "RE: %lf\n", RE);
  // fprintf(stderr, "CA: %lf\n", CA);
  // fprintf(stderr, "LX: %lf\n", LX);
  // fprintf(stderr, "TH: %lf\n", THETA);
  // fprintf(stderr, " N: %d\n", N);
  // fprintf(stderr, " M: %d\n", M);
  // fprintf(stderr, " P: %d\n", P);

  /* compute J, Psi, Phi, and K0 */
  /* Jacobian */
  double **J = malloc_f2d(2*N, 2*N);
  wr_jacobian(J);
  // fprintf(stderr, "J\n  %10lf %10lf\n  %10lf %10lf\n", J[0+N][0+N], J[0+N][1+N], J[1+N][0+N], J[1+N][1+N]);

  /* actuator matrix */
  double **Psi = malloc_f2d(2*N, M);
  wr_actuator(Psi);
  // fprintf(stderr, "Psi\n  %10lf %10lf\n  %10lf %10lf\n", Psi[0][0], Psi[0][1], Psi[1][0], Psi[1][1]);

  /* observer matrix (actually the transpose) */
  double **Phi = malloc_f2d(2*N, P);
  wr_observer(Phi);
  // fprintf(stderr, "Phi\n  %10lf %10lf\n  %10lf %10lf\n", Phi[0][0], Phi[0][1], Phi[1][0], Phi[1][1]);

  /* ======================================= */
  /* find an initial stabilising guess for K */
  /* ======================================= */
  // following Ilka 2023
  // const double vinv =  1.0/(1.0-MU);

  // /* compute Stil = Psi * Vinv * Psi' */
  // double **Stil = malloc_f2d(2*N, 2*N);
  // for (int i = 0; i < 2*N; i++) {
  //   for (int j = 0; j < 2*N; j++) {
  //     Stil[i][j] = 0.0;
  //     for (int k = 0; k < M; k++) {
  //       Stil[i][j] += Psi[i][k] * Psi[j][k] * vinv;
  //     } // k end
  //   } // j end
  // } // i end

  // // fprintf(stderr, "Stil\n  %10lf %10lf\n  %10lf %10lf\n", Stil[0][0], Stil[0][1], Stil[1][0], Stil[1][1]);

  // /* compute CC = Cp*C - I */
  // double **CC = malloc_f2d(2*N, 2*N);
  // for (int i = 0; i < 2*N; i++) {
  //   for (int j = 0; j < 2*N; j++) {
  //     CC[i][j] = 0.0;
  //   } // j end
  //   CC[i][i] -= 1.0;
  // } // i end
  // for (int i = 0; i < N; i++) {
  //   for (int j = 0; j < N; j++) {
  //     for (int k = 0; k < P; k++) {
  //       CC[i][j] += Phi[i][k] * Phi[j][k];
  //     } // k end
  //   } // j end
  // } // i end

  // // fprintf(stderr, "CC\n  %10lf %10lf\n  %10lf %10lf\n", CC[0][0], CC[0][1], CC[1][0], CC[1][1]);

  // /* compute initial guess for P */
  // // TODO: might need to multiply by -1
  // double **Pj = malloc_f2d(2*N, 2*N);
  // dlqr(J, Stil, 1.0, 0.1, 2*N, 2*N, Pj);

  // // for (int i = 0; i < N; i++) {
  // //   for (int j = 0; j < N; j++) {
  // //     Pj[i][j] *= -1.0;
  // //   } // j end
  // // } // i end

  // // fprintf(stderr, "Pj\n  %10lf %10lf\n  %10lf %10lf\n", Pj[0][0], Pj[0][1], Pj[1][0], Pj[1][1]);

  // // recompute J and Stil
  // // TODO: required?
  // wr_jacobian(J);
  // for (int i = 0; i < 2*N; i++) {
  //   for (int j = 0; j < 2*N; j++) {
  //     Stil[i][j] = 0.0;
  //     for (int k = 0; k < M; k++) {
  //       Stil[i][j] += Psi[i][k] * Psi[j][k] * vinv;
  //     } // k end
  //   } // j end
  // } // i end


  // /* iterate */
  // const int iternum = 10;
  // double **Gj = malloc_f2d(2*N, 2*N);
  // double **RPj = malloc_f2d(2*N, 2*N);
  // double **BPj = malloc_f2d(M, 2*N);
  // double **PjB = malloc_f2d(2*N, M);
  // for (int l = 0; l < iternum; l++) {
  //   /* Gj = Rinv*(B'*Pj)*CC */
  //   for (int i = 0; i < M; i++) {
  //     for (int j = 0; j < 2*N; j++) {
  //       Gj[i][j] = 0.0;
  //       for (int k = 0; k < 2*N; k++) {
  //         for (int m = 0; m < 2*N; m++) {
  //           Gj[i][j] += Psi[m][i] * Pj[m][k] * CC[k][j] * vinv;
  //         } // m end
  //       } // k end
  //     } // j end
  //   } // i end

  //   // fprintf(stderr, "Gj\n  %10lf %10lf\n  %10lf %10lf\n", Gj[0][0], Gj[0][1], Gj[1][0], Gj[1][1]);

  //   /* B' * Pj */
  //   for (int i = 0; i < M; i++) {
  //     for (int j = 0; j < 2*N; j++) {
  //       BPj[i][j] = 0.0;
  //       for (int k = 0; k < 2*N; k++) {
  //         BPj[i][j] += Psi[k][i] * Pj[k][j];
  //       } // k end
  //     } // j end
  //   } // i end

  //   /* B' * Pj */
  //   for (int i = 0; i < 2*N; i++) {
  //     for (int j = 0; j < M; j++) {
  //       PjB[i][j] = 0.0;
  //       for (int k = 0; k < 2*N; k++) {
  //         PjB[i][j] += Pj[i][k] * Psi[k][j];
  //       } // k end
  //     } // j end
  //   } // i end

  //   /* RPj = Gj'*R*Gj + A'*Pj + Pj*A + Q - Rinv*(Pj*B)*(B'*Pj) */
  //   for (int i = 0; i < 2*N; i++) {
  //     for (int j = 0; j < 2*N; j++) {
  //       RPj[i][j] = 0.0;
  //       for (int k = 0; k < 2*N; k++) {
  //         RPj[i][j] += Gj[k][i]*Gj[k][j]*(1.0-MU) + J[k][i]*Pj[k][j] + Pj[i][k]*J[k][j];
  //       } // k end
  //       for (int k = 0; k < M; k++) {
  //         RPj[i][j] -= vinv*PjB[i][k]*PjB[k][j];
  //       } // k end
  //     } // j end
  //   } // i end

  //   // fprintf(stderr, "RPj\n  %10lf %10lf\n  %10lf %10lf\n", RPj[0][0], RPj[0][1], RPj[1][0], RPj[1][1]);

  //   /* compute (Frobenius) norm of RPj */
  //   // TODO: maybe don't bother!
  //   double norm = LAPACKE_dlange(LAPACK_ROW_MAJOR, 'F', 2*N, 2*N, *RPj, 2*N);
  //   fprintf(stderr, "tr: %lf\n", norm);
  //   if (norm < 1.0e-6) {
  //     break;
  //   } else {
  //     /* compute Ac (reusing space from Gj) */
  //     for (int i = 0; i < 2*N; i++) {
  //       for (int j = 0; j < 2*N; j++) {
  //         Gj[i][j] = J[i][j];
  //         for (int k = 0; k < 2*N; k++) {
  //           Gj[i][j] -= Stil[i][k] * Pj[k][j];
  //         } // k end
  //       } // j end
  //     } // i end

  //     // fprintf(stderr, "Ac\n  %10lf %10lf\n  %10lf %10lf\n", Gj[0][0], Gj[0][1], Gj[1][0], Gj[1][1]);

  //     /* solve for step (Xj = RPj after dlyap called) */
  //     dlyap(Gj, RPj, 2*N, 1);

  //     /* step */
  //     for (int i = 0; i < 2*N; i++) {
  //       for (int j = 0; j < 2*N; j++) {
  //         Pj[i][j] += RPj[i][j];
  //       } // j end
  //     } // i end

  //     // fprintf(stderr, "RPj\n  %10lf %10lf\n  %10lf %10lf\n", RPj[0][0], RPj[0][1], RPj[1][0], RPj[1][1]);
  //     // fprintf(stderr, "Gj\n  %10lf %10lf\n  %10lf %10lf\n", Gj[0][0], Gj[0][1], Gj[1][0], Gj[1][1]);
  //     // ABORT("expected end");
  //   }
  // } // l end

  // /* compute initial guess */
  // /* F = Rinv*(B'*Pj)*Cp */
  // for (int i = 0; i < M; i++) {
  //   for (int j = 0; j < P; j++) {
  //     STATIC_K[i][j] = 0.0;
  //     for (int k = 0; k < 2*N; k++) {
  //       for (int m = 0; m < 2*N; m++) {
  //         STATIC_K[i][j] += Psi[m][i] * Pj[m][k] * Phi[k][j] * vinv;
  //       } // m end
  //     } // k end
  //   } // j end
  // } // i end

  // free_2d(Stil);
  // free_2d(CC);
  // free_2d(Pj);
  // free_2d(Gj);
  // free_2d(RPj);

  FILE *fp = fopen("k0.dat", "r");
  if (!fp) { ABORT("failed to open k0.dat"); }
  for (int i = 0; i < M; i++) {
    for (int j = 0; j < P; j++) {
      fscanf(fp, "%lf", &(STATIC_K[i][j]));
    } // j end
  } // i end
  // fread(*STATIC_K, sizeof(double), M*P, fp);
  fclose(fp);

  fprintf(stderr, "READ SUCCESSFUL\n");
  fprintf(stderr, "STATIC_K\n  %15lf %15lf\n  %15lf %15lf\n", STATIC_K[0][0], STATIC_K[0][1], STATIC_K[1][0], STATIC_K[1][1]);
  fflush(stderr);

  // /* set up matrices for iteration */
  // double **Ak = malloc_f2d(2*N, 2*N);
  double **KPHIk = malloc_f2d(M, 2*N);


  // /* resulting K*Phi */
  // for (int i = 0; i < M; i++) {
  //   for (int j = 0; j < 2*N; j++) {
  //     KPHIk[i][j] = 0.0;
  //     for (int k = 0; k < P; k++) {
  //       KPHIk[i][j] += STATIC_K[i][k] * Phi[j][k];
  //     } // k end
  //   } // j end
  // } // i end

  // /* closed-loop system matrix */
  // for (int i = 0; i < 2*N; i++) {
  //   for (int j = 0; j < 2*N; j++) {
  //     Ak[i][j] = J[i][j];
  //     for (int k = 0; k < M; k++) {
  //       Ak[i][j] -= Psi[i][k] * KPHIk[k][j];
  //     } // k end
  //   } // j end
  // } // i end

  // fprintf(stderr, "Ak\n  %15lf %15lf\n  %15lf %15lf\n", Ak[N+0][N+0], Ak[N+0][N+1], Ak[N+1][N+0], Ak[N+1][N+1]);


  // /* check that K0 is stable */
  // double lmax;
  // dsr(Ak, N, &lmax);
  // if (lmax >= 0.0) {
  //   // ABORT("K0 is not stable (l = %lf)", lmax);
  //   fprintf(stderr, "K0 is not stable (l = %lf)", lmax);
  // } else {
  //   fprintf(stderr, "K0 is stable (l = %lf)", lmax);
  // }

  // /* iterate to find K */
  // double c = 0.0; // cost
  // double dc; // cost difference
  // double iv;
  // double **Pk = malloc_f2d(2*N, 2*N);
  // double **Sk = malloc_f2d(2*N, 2*N);

  // double **K1 = malloc_f2d(M, P); // search direction
  // double **W1 = malloc_f2d(M, 2*N); // TODO: remove if possible
  // double **W2 = malloc_f2d(2*N, P); // TODO: remove if possible
  // double **W3 = malloc_f2d(P, P); // TODO: remove if possible
  // int *ipiv = malloc(P*sizeof(int));
  // while (1) {
  //   /* K * Phi */
  //   for (int i = 0; i < M; i++) {
  //     for (int j = 0; j < 2*N; j++) {
  //       KPHIk[i][j] = 0.0;
  //       for (int k = 0; k < P; k++) {
  //         KPHIk[i][j] += STATIC_K[i][k] * Phi[j][k];
  //       } // k end
  //     } // j end
  //   } // i end


  //   /* Ak = (J - Psi * K * Phi)' */
  //   for (int i = 0; i < 2*N; i++) {
  //     for (int j = 0; j < 2*N; j++) {
  //       Ak[j][i] = J[i][j];
  //       for (int k = 0; k < M; k++) {
  //         Ak[j][i] -= Psi[i][k] * KPHIk[k][j];
  //       } // k end
  //     } // j end
  //   } // i end

  //   /* Q = v * (K * Phi)' * K * Phi + u * I */
  //   for (int i = 0; i < 2*N; i++) { //
  //     for (int j = 0; j < 2*N; j++) {
  //       Pk[i][j] = 0.0;
  //       for (int k = 0; k < M; k++) {
  //         Pk[i][j] += (1-MU) * KPHIk[k][i] * KPHIk[k][j];
  //       } // k end
  //     } // j end
  //     Pk[i][i] += DX * MU; // TODO: check the factor of DX
  //   } // i end

  //   /* solve for Pk */
  //   dlyap(Ak, Pk, 2*N, 0);


  //   /* Ak = J - Psi * K * Phi */
  //   for (int i = 0; i < 2*N; i++) {
  //     for (int j = 0; j < 2*N; j++) {
  //       Ak[i][j] = J[i][j];
  //       for (int k = 0; k < M; k++) {
  //         Ak[i][j] -= Psi[i][k] * KPHIk[k][j];
  //       } // k end
  //     } // j end
  //   } // i end

  //   /* Q = I */
  //   for (int i = 0; i < 2*N; i++) {
  //     for (int j = 0; j < 2*N; j++) {
  //       Sk[i][j] = 0.0;
  //     } // j end
  //     Sk[i][i] = 1.0;
  //   } // i end

  //   /* solve for Sk */
  //   // TODO: this can be sped up singe Sk = I
  //   dlyap(Ak, Sk, 2*N, 0);


  //   /* compute delta cost */
  //   dc = 0.5 * dtr(Pk, 2*N) - c;
  //   c = c + dc;
  //   fprintf(stderr, "dc: %lf\n", dc);


  //   /* compute gain update direction */
  //   iv = 1/(1-MU); // scale factor
  //   for (int i = 0; i < M; i++) { // interim W1
  //     for (int j = 0; j < 2*N; j++) {
  //       W1[i][j] = 0.0;
  //       for (int k = 0; k < 2*N; k++) {
  //         W1[i][j] += iv * Psi[k][i] * Pk[k][j];
  //       } // k end
  //     } // j end
  //   } // i end
  //   for (int i = 0; i < 2*N; i++) { // interim W2
  //     for (int j = 0; j < P; j++) {
  //       W2[i][j] = 0.0;
  //       for (int k = 0; k < 2*N; k++) {
  //         W2[i][j] += Sk[i][k] * Phi[k][j];
  //       } // k end
  //     } // j end
  //   } // i end
  //   for (int i = 0; i < M; i++) { // pre-solve K = W1 * W2
  //     for (int j = 0; j < P; j++) {
  //       K1[i][j] = 0.0;
  //       for (int k = 0; k < 2*N; k++) {
  //         K1[i][j] += W1[i][k] * W2[k][j];
  //       } // k end
  //     } // j end
  //   } // i end
  //   for (int i = 0; i < P; i++) { // interim W3 (to invert)
  //     for (int j = 0; j < P; j++) {
  //       W3[i][j] = 0.0;
  //       for (int k = 0; k < 2*N; k++) {
  //         W3[i][j] += Phi[k][i] * W2[k][j];
  //       } // k end
  //     } // j end
  //   } // i end

  //   // solve W3 * K1 = K1 (ie the transpose of the problem)
  //   LAPACKE_dgesv(LAPACK_COL_MAJOR, P, M, *W3, P, ipiv, *K1, P);


  //   /* update K */
  //   const double a = 0.005; // TODO: set optimally
  //   for (int i = 0; i < M; i++) { // pre-solve K = W1 * W2
  //     for (int j = 0; j < P; j++) {
  //       STATIC_K[i][j] = (1-a)*STATIC_K[i][j] + a*K1[i][j]; // note the transpose
  //     } // j end
  //   } // i end


  //   /* finish if dC is small enough */ // TODO: how small is enough
  //   if (fabs(dc) < 1e-4) {
  //     break;
  //   }
  // }

  /* K * Phi */
  for (int i = 0; i < M; i++) {
    for (int j = 0; j < 2*N; j++) {
      KPHIk[i][j] = 0.0;
      for (int k = 0; k < P; k++) {
        KPHIk[i][j] += STATIC_K[i][k] * Phi[j][k];
      } // k end
    } // j end
  } // i end

  /* apply flux approximation */
  for (int i = 0; i < M; i++) {
    for (int j = 0; j < N; j++) {
      STATIC_KPHI[i][j] = KPHIk[i][j] + 2/3.0 * KPHIk[i][N+j];
    } // j end
  } // i end


  free_2d(J);
  free_2d(Psi);
  free_2d(Phi);
  // free_2d(Ak);
  // free_2d(Pk);
  // free_2d(Sk);
  // free_2d(K1);
  // free_2d(W1);
  // free_2d(W2);
  // free_2d(W3);
  // free(ipiv);
}


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
      static_wr_compute_KPHI(STATIC_KPHI);
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
