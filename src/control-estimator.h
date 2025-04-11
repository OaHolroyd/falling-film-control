#ifndef CONTROL_ESTIMATOR_H
#define CONTROL_ESTIMATOR_H

#include <math.h>

#include "c-utils.h"
#include "control-core.h"
#include "linalg.h"
#include "wr.h"

#include <string.h>

static double *EST_y;   // difference in observation between real and estimator
static double *EST_f;   // forcing term
static double *EST_ff;  // forcing term (of the main system)
static double *EST_h;   // height estimate
static double *EST_q;   // flux estimate
static double *EST_h0;  // height estimate (previous time step)
static double *EST_q0;  // flux estimate (previous time step)
static double *EST_res; // residual (for time stepping)
static double **EST_L;  // estimator forcing matrix
static double **EST_K;  // main forcing matrix
static double **EST_C;  // observer matrix

// storage used in inverting the Jacobian
static double *EST_cl2;
static double *EST_cl1;
static double *EST_cd0;
static double *EST_cu1;
static double *EST_sl2;
static double *EST_sl1;
static double *EST_sd0;
static double *EST_su1;
static double *EST_su2;
static double *EST_k0;
static double *EST_k1;
static double *EST_work_z;

static double *EST_work;

// derivative macros
#define WRAP(i) ((i + N) % N)

// centre-to-centre spacing
#define D1C(z, i) (0.5 * ((z[WRAP(i + 1)] - z[WRAP(i - 1)]) / DX))

// left and right spacing
#define D0L(z, i) (0.5 * (z[WRAP(i + 1)] + z[WRAP(i)]))
#define D1L(z, i) ((z[WRAP(i + 1)] - z[WRAP(i)]) / DX)
#define D0R(z, i) (0.5 * (z[WRAP(i)] + z[WRAP(i - 1)]))
#define D1R(z, i) ((z[WRAP(i)] - z[WRAP(i - 1)]) / DX)
#define D3R(z, i)                                                              \
  ((-z[WRAP(i - 2)] + 3.0 * z[WRAP(i - 1)] - 3.0 * z[WRAP(i)] +                \
    z[WRAP(i + 1)]) /                                                          \
   (DX * DX * DX))

/* ========================================================================== */
/*   AUXILIARY FUNCTION DEFINITIONS                                           */
/* ========================================================================== */
/* solve the (transpose) LQR problem to compute L. This requires EST_C to have
 * been filled with the Benney or WR observer. */
void est_forcing_matrix(double **L) {
  // transpose of L
  double **Lt = malloc_f2d(P, 2 * N);

  // transpose of the system matrix
  double **At = malloc_f2d(2 * N, 2 * N);
  wr_jacobian(At);
  for (int i = 0; i < 2 * N; i++) {
    for (int j = 0; j < 2 * N; j++) {
      double tmp = At[i][j];
      At[i][j] = At[j][i];
      At[j][i] = tmp;
    } // j end
  } // i end

  // transpose of the observer matrix
  // TODO: maybe have this be self contained
  double **Ct = malloc_f2d(2 * N, P);
  // TODO use this to zero out the second half
  // memset(Ct[N], 0, N * P * sizeof(double));
  for (int i = 0; i < 2 * N; i++) {
    for (int j = 0; j < P; j++) {
      Ct[i][j] = 0.0;
    } // j end
  }
  for (int i = 0; i < N; i++) { // only allowed to observe the height
    for (int j = 0; j < P; j++) {
      Ct[i][j] = EST_C[i][j];
    } // j end
  }

  /* compute L using LQR */
  double u = 1.0 / sqrt(DX); // cost of estimator error
  double v = 1.0 * sqrt(DX); // cost of control (should scale with DX)
  dlqr(At, Ct, u, v, 2 * N, P, Lt);

  // transpose Lt to get L
  for (int i = 0; i < 2 * N; i++) {
    for (int j = 0; j < P; j++) {
      L[i][j] = Lt[j][i];
    } // j end
  }

  free_2d(Lt);
  free_2d(At);
  free_2d(Ct);
}

/* compute the gain matrix mapping the estimator to actuator strengths for the
 * main problem */
void est_gain_matrix(double **K) {
  /* Jacobian */
  double **A = malloc_f2d(2 * N, 2 * N);
  wr_jacobian(A);

  /* actuator matrix */
  double **B = malloc_f2d(2 * N, M);
  wr_actuator(B);

  /* full control matrix */
  dlqr(A, B, sqrt(DX) * MU, 1.0 / sqrt(DX), 2 * N, M, K);

  free_2d(A);
  free_2d(B);
}

/* step the estimator forward in time using observations of the real height H.
 * Returns the number of iterations used. */
int est_update(double dt, double *H) {
  // don't bother doing anything if EST_h < 0.0 or is nan
  for (int i = 0; i < N; i++) {
    if (EST_h[i] < 0.0 || isnan(EST_h[i])) {
      return 0;
    }
  }

  // work out the difference in observed height and the estimate
  for (int i = 0; i < P; i++) {
    EST_y[i] = 0.0;
    for (int j = 0; j < N; j++) {
      EST_y[i] += EST_C[j][i] * (H[j] - EST_h[j]);
    } // j end
  }

  // compute forcing term and rhs
  for (int i = 0; i < 2 * N; i++) {
    EST_f[i] = 0.0;
    for (int j = 0; j < P; j++) {
      EST_f[i] += EST_L[i][j] * EST_y[j];
    } // j end
  }

  // compute the forcing term for the main system
  for (int i = 0; i < N; i++) {
    EST_ff[i] = control(ITOX(i));
  }

  // pack the data into the wrapper struct
  struct wr_data data = {
      .n = N,
      .dx = DX,
      .theta = THETA,
      .re = RE,
      .ca = CA,
      .h = EST_h,
      .q = EST_q,
      .h0 = EST_h0,
      .q0 = EST_q0,
      .fa = EST_ff,
      .f = EST_f,
      .work = EST_work,
  };

  // implicit time-stepping (for stability)
  // iterate to the solution for the next timestep
  return wr_step(&data, dt, 1.0e-6, 100);
}

/* ========================================================================== */
/*   FUNCTION DEFINITIONS                                                     */
/* ========================================================================== */
/* [REQUIRED] internal setup */
void est_set(void) {
  // height observations
  EST_y = malloc(P * sizeof(double));

  // (transpose of the) observer matrix
  EST_C = malloc_f2d(N, P);
  benney_observer(EST_C);

  // estimator forcing matrix
  EST_L = malloc_f2d(2 * N, P);
  est_forcing_matrix(EST_L);

  // main forcing matrix
  EST_K = malloc_f2d(M, 2 * N);
  est_gain_matrix(EST_K);

  // forcing term
  EST_f = malloc(2 * N * sizeof(double));
  EST_ff = malloc(N * sizeof(double));

  // height estimator
  EST_h = malloc(N * sizeof(double));
  EST_h0 = malloc(N * sizeof(double));
  for (int i = 0; i < N; i++) {
    EST_h[i] = 1.0;
    EST_h0[i] = 1.0;
  }

  // flux estimator
  EST_q = malloc(N * sizeof(double));
  EST_q0 = malloc(N * sizeof(double));
  for (int i = 0; i < N; i++) {
    EST_q[i] = 2.0 / 3.0;
    EST_q0[i] = 2.0 / 3.0;
  }

  // residual and Jacobian (for time stepping)
  EST_res = malloc(2 * N * sizeof(double));

  EST_cl2 = malloc(N * sizeof(double));
  EST_cl1 = malloc(N * sizeof(double));
  EST_cd0 = malloc(N * sizeof(double));
  EST_cu1 = malloc(N * sizeof(double));
  EST_sl2 = malloc(N * sizeof(double));
  EST_sl1 = malloc(N * sizeof(double));
  EST_sd0 = malloc(N * sizeof(double));
  EST_su1 = malloc(N * sizeof(double));
  EST_su2 = malloc(N * sizeof(double));
  EST_k0 = malloc(N * sizeof(double));
  EST_k1 = malloc(N * sizeof(double));
  EST_work_z = malloc(N * sizeof(double));

  EST_work = malloc(14 * N * sizeof(double));

  /* pick from the available ROMs */
  switch (RT) {
  // case BENNEY:
  //      lqr_benney_compute_K(LQR_K);
  // break;
  case WR:
    //      lqr_wr_compute_K(LQR_K);
    break;
  default:
    ABORT("invalid ROM type %d", RT);
  }
}

/* [REQUIRED] internal free */
void est_free(void) {
  free(EST_h);
  free(EST_h0);
  free(EST_q);
  free(EST_q0);
  free(EST_res);
  free(EST_f);
  free(EST_ff);
  free(EST_y);
  free_2d(EST_L);
  free_2d(EST_K);
  free_2d(EST_C);

  free(EST_cl2);
  free(EST_cl1);
  free(EST_cd0);
  free(EST_cu1);
  free(EST_sl2);
  free(EST_sl1);
  free(EST_sd0);
  free(EST_su1);
  free(EST_su2);
  free(EST_k0);
  free(EST_k1);
  free(EST_work_z);
}

/* [REQUIRED] steps the system forward in time given the interfacial height */
int est_step(double dt, double *h, double *q, int control_on) {
  /* u = K * (h-1) */
  for (int i = 0; i < M; i++) {
    Amag[i] = 0.0;

    if (control_on) {
      for (int j = 0; j < N; j++) {
        Amag[i] += EST_K[i][j] * (EST_h[j] - 1.0);
        Amag[i] += EST_K[i][j + N] * (D0L(EST_q, j) - 2.0 / 3.0);
      } // j end
    }
  } // i end

  /* update the estimator */
  return est_update(dt, h);
}

/* [REQUIRED] returns the estimator as a function of x */
double est_estimator(double x) { return interp(x, EST_h); }

/* [REQUIRED] outputs the internal matrices */
void est_output(void) {
  // char fname[128];
  // sprintf(fname, "out/L_%d.dat", N);
  // output_d2d(fname, EST_L, 2 * N, P);
  // sprintf(fname, "out/K_%d.dat", N);
  // output_d2d(fname, EST_K, M, 2 * N);

  output_d2d("out/L.dat", EST_L, 2 * N, P);
  output_d2d("out/K.dat", EST_K, M, 2 * N);
}

/* [REQUIRED] generates the control matrix CM = F*K */
void est_matrix(double **CM) {
  /* forcing matrix */
  double **F = malloc_f2d(N, M);
  forcing_matrix(F);

  for (int i = 0; i < N; i++) {
    for (int j = 0; j < N; j++) {
      CM[i][j] = 0.0;
      for (int k = 0; k < M; k++) {
        //        CM[i][j] += F[i][k]*LQR_K[k][j];
      } // k end
    } // j end
  } // i end

  free_2d(F);
}

#endif
