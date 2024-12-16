#ifndef CONTROL_QQR_H
#define CONTROL_QQR_H

#include <math.h>

#include "c-utils.h"
#include "linalg.h"
#include "control-core.h"


static double **QQR_C0; /* control operator */
static double **QQR_C1; /* control operator */


/* ========================================================================== */
/*   AUXILIARY FUNCTION DEFINITIONS                                           */
/* ========================================================================== */
/* compute the control matrix in the Benney case */
void qqr_benney_compute_matrices(double **qqr_c0, double **qqr_c1) {
  /* Jacobian */
  double **J = malloc_f2d(N, N);
  benney_jacobian(J);

  /* actuator matrix */
  double **Psi = malloc_f2d(N, M);
  benney_actuator(Psi);

  /* control matrix */
  dlqr(J, Psi, DX*MU, 1-MU, N, M, qqr_c0);

  free_2d(J);
  free_2d(Psi);
}

/* compute the control matrix in the weighted-residuals case */
void qqr_wr_compute_matrices(double **qqr_c0, double **qqr_c1) {
  /* Jacobian */
  double **A = malloc_f2d(2*N, 2*N);
  wr_jacobian(A);
  output_d2d("out/A.dat", A, 2*N, 2*N);

  /* actuator matrix */
  double **B = malloc_f2d(2*N, M);
  wr_actuator(B);
  output_d2d("out/B.dat", B, 2*N, M);

  /* linear control matrix */
  dlqr(A, B, DX*MU, 1-MU, 2*N, M, qqr_c0);

  /* work out quadratic control matrix */
  double complex **QQR_P = malloc_z2d(2*N, 2*N);
  riccati(A, B, DX*MU, 1-MU, 2*N, M, QQR_P[0]);
  output_z2d("out/P.dat", QQR_P, 2*N, 2*N);

  // BBt = B * B'
  double **BBt = malloc_f2d(2*N, 2*N);
  for (int i = 0; i < 2*N; i++) {
    for (int j = 0; j < 2*N; j++) {
      BBt[i][j] = 0.0;
      for (int k = 0; k < M; k++) {
        BBt[i][j] += B[i][k] * B[j][k];
      }
    }
  }
  output_d2d("out/BBt.dat", BBt, 2*N, 2*N);

  // C = A' - 1/(1-MU) * P*B*B';
  double **C = malloc_f2d(2*N, 2*N);
  for (int i = 0; i < 2*N; i++) {
    for (int j = 0; j < 2*N; j++) {
      C[i][j] = 0.0;
      for (int k = 0; k < 2*N; k++) {
        C[i][j] += creal(QQR_P[i][k]) * BBt[k][j];
      }
      C[i][j] *= -1.0 / (1.0 - MU);
      C[i][j] += A[j][i];
    }
  }
  output_d2d("out/C.dat", C, 2*N, 2*N);

  // D = C \ P = [A' - 1/(1-MU) * P*B*B']^-1 * P
  double **D = malloc_f2d(2*N, 2*N);
  for (int i = 0; i < 2*N; i++) {
    for (int j = 0; j < 2*N; j++) {
      D[i][j] = creal(QQR_P[i][j]);
    }
  }
  output_d2d("out/Preal.dat", D, 2*N, 2*N);

  dgesv(C, D, 2*N);
  output_d2d("out/D.dat", D, 2*N, 2*N);

  // qqr_c1 = -1/(1-MU) * B' * D = -1/(1-MU) * B' * [A' - 1/(1-MU) * P*B*B'] \ P
  for (int i = 0; i < M; i++) {
    for (int j = 0; j < 2*N; j++) {
      qqr_c1[i][j] = 0.0;
      for (int k = 0; k < 2*N; k++) {
        qqr_c1[i][j] += B[k][i] * D[k][j];
      } // k end
      qqr_c1[i][j] *= -1.0 / (1.0 - MU);
    } // j end
  } // i end
  output_d2d("out/C1.dat", qqr_c1, M, 2*N);
  // exit(1);


  free_2d(A);
  free_2d(B);
  free_2d(QQR_P);
  free_2d(BBt);
  free_2d(C);
  free_2d(D);
}


/* ========================================================================== */
/*   FUNCTION DEFINITIONS                                                     */
/* ========================================================================== */
/* [REQUIRED] internal setup */
void qqr_set(void) {
  QQR_C0 = malloc_f2d(M, 2*N);
  QQR_C1 = malloc_f2d(M, 2*N);

  /* pick from the available ROMs */
  switch (RT) {
    case BENNEY:
      qqr_benney_compute_matrices(QQR_C0, QQR_C1);
      break;
    case WR:
      qqr_wr_compute_matrices(QQR_C0, QQR_C1);
      break;
    default :
      ABORT("invalid ROM type %d", RT);
  }
}

/* [REQUIRED] internal free */
void qqr_free(void) {
  free(QQR_C0);
  free(QQR_C1);
}

/* [REQUIRED] steps the system forward in time given the interfacial height */
void qqr_step(double dt, double *h, double *q) {
  double H[5];
  double Q[3];

  /* f = K * (h-1) */
  for (int i = 0; i < M; i++) {
    Amag[i] = 0.0;

    /* height component */
    for (int j = 0; j < N; j++) {
      const double hj = h[j] - 1.0;
      const double hhj = 0.0;

      Amag[i] += QQR_C0[i][j] * hj + QQR_C1[i][j] * hhj;
    } // j end

    /* flux component */
    for (int j = 0; j < N; j++) {
      /* get data at j and its neighbours */
      if (j == 0) {
        H[0] = h[N-2] - 1.0;
        H[1] = h[N-1] - 1.0;
        H[2] = h[0] - 1.0;
        H[3] = h[1] - 1.0;
        H[4] = h[2] - 1.0;

        Q[0] = q[N-1] - 2.0/3.0;
        Q[1] = q[0] - 2.0/3.0;
        Q[2] = q[1] - 2.0/3.0;
      } else if (j == 1) {
        H[0] = h[N-1] - 1.0;
        H[1] = h[0] - 1.0;
        H[2] = h[1] - 1.0;
        H[3] = h[2] - 1.0;
        H[4] = h[3] - 1.0;

        Q[0] = q[0] - 2.0/3.0;
        Q[1] = q[1] - 2.0/3.0;
        Q[2] = q[2] - 2.0/3.0;
      } else if (j == N-2) {
        H[0] = h[N-4] - 1.0;
        H[1] = h[N-3] - 1.0;
        H[2] = h[N-2] - 1.0;
        H[3] = h[N-1] - 1.0;
        H[4] = h[0] - 1.0;

        Q[0] = q[N-3] - 2.0/3.0;
        Q[1] = q[N-2] - 2.0/3.0;
        Q[2] = q[N-1] - 2.0/3.0;
      } else if (j == N-1) {
        H[0] = h[N-3] - 1.0;
        H[1] = h[N-2] - 1.0;
        H[2] = h[N-1] - 1.0;
        H[3] = h[0] - 1.0;
        H[4] = h[1] - 1.0;

        Q[0] = q[N-2] - 2.0/3.0;
        Q[1] = q[N-1] - 2.0/3.0;
        Q[2] = q[0] - 2.0/3.0;
      } else {
        H[0] = h[j-2] - 1.0;
        H[1] = h[j-1] - 1.0;
        H[2] = h[j] - 1.0;
        H[3] = h[j+1] - 1.0;
        H[4] = h[j+2] - 1.0;

        Q[0] = q[j-1] - 2.0/3.0;
        Q[1] = q[j] - 2.0/3.0;
        Q[2] = q[j+1] - 2.0/3.0;
      }

      /* compute derivatives */
      const double hx = 0.5 * (H[3]-H[1]) / DX;
      const double hxxx = 0.5 * (H[4] - 2*H[3] + 2*H[1] - H[0]) / (DX*DX*DX);
      const double qx = 0.5 * (Q[2]-Q[0]) / DX;

      // TODO: compute nonlinear bit correctly
      const double hj = h[j] - 1.0;
      const double qj = q[j] - 2.0/3.0;
      const double qqj = -(72.0*RE*qj*hx-68.0*RE*hj*qx-102.0*RE*qj*qx+120*hj*hj-210*hj*hx+105*hj*hxxx/CA)/105.0;

      Amag[i] += QQR_C0[i][N+j] * qj + DEL * QQR_C1[i][N+j] * qqj;
    } // j end
  } // i end
}

/* [REQUIRED] returns the estimator as a function of x */
double qqr_estimator(double x) {

  return 0.0;
}

/* [REQUIRED] outputs the internal matrices */
void qqr_output(void) {
  output_d2d("out/C0.dat", QQR_C0, M, 2*N);
  output_d2d("out/C1.dat", QQR_C1, M, 2*N);
}

/* [REQUIRED] generates the control matrix CM = F*K */
void qqr_matrix(double **CM) {
  exit(1);

  /* forcing matrix */
  double **F = malloc_f2d(N, M);
  forcing_matrix(F);

  for (int i = 0; i < N; i++) {
    for (int j = 0; j < N; j++) {
      CM[i][j] = 0.0;
      for (int k = 0; k < M; k++) {
        CM[i][j] += F[i][k]*QQR_C0[k][j];
      } // k end
    } // j end
  } // i end

  free_2d(F);
}

#endif
