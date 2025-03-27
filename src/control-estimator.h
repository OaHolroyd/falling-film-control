#ifndef CONTROL_ESTIMATOR_H
#define CONTROL_ESTIMATOR_H

#include <math.h>

#include "c-utils.h"
#include "control-core.h"
#include "linalg.h"

static double *EST_y;   // difference in observation between real and estimator
static double *EST_f;   // forcing term
static double *EST_ff;  // forcing term (of the main system)
static double *EST_h;   // height estimate
static double *EST_q;   // flux estimate
static double *EST_h0;  // height estimate (previous time step)
static double *EST_q0;  // flux estimate (previous time step)
static double *EST_res; // residual (for time stepping)
static double **EST_J;  // Jacobian (for time stepping)
static double **EST_L;  // estimator forcing matrix
static double **EST_K;  // main forcing matrix
static double **EST_C;  // observer matrix

// derivative macros
#define WRAP(i) ((i + N) % N)
#define D1(z, i) (0.5 * ((z[WRAP(i + 1)] - z[WRAP(i - 1)]) / DX))
#define D2(z, i)                                                               \
  ((z[WRAP(i + 1)] - 2.0 * z[WRAP(i)] + z[WRAP(i - 1)]) / (DX * DX))
#define D3(z, i)                                                               \
  ((z[WRAP(i + 2)] - 2.0 * z[WRAP(i + 1)] + 2.0 * z[WRAP(i - 1)] -             \
    z[WRAP(i - 2)]) /                                                          \
   (2.0 * DX * DX * DX))
#define D4(z, i)                                                               \
  ((z[WRAP(i + 2)] - 4.0 * z[WRAP(i + 1)] + 6.0 * z[WRAP(i)] -                 \
    4.0 * z[WRAP(i - 1)] + z[WRAP(i - 2)]) /                                   \
   (DX * DX * DX * DX))

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
  double **Ct = malloc_f2d(2 * N, P);
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
  double u = 1.0; // cost of estimator error
  double v = 0.1; // cost of control effort (don't care)
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

void est_gain_matrix(double **K) {
  /* Jacobian */
  double **A = malloc_f2d(2 * N, 2 * N);
  wr_jacobian(A);

  /* actuator matrix */
  double **B = malloc_f2d(2 * N, M);
  wr_actuator(B);

  /* full control matrix */
  dlqr(A, B, DX * MU, 1 - MU, 2 * N, M, K);

  free_2d(A);
  free_2d(B);
}

/* compute the residual and return the square of the norm */
double est_compute_residual(double dt, double *res) {
  double res_norm_2 = 0.0;
  double *h = EST_h;
  double *h0 = EST_h0;
  double *q = EST_q;
  double *q0 = EST_q0;
  double *f = EST_f;
  double *ff = EST_ff;

  for (int i = 0; i < N; i++) {
    // H component
    res[i] = 2.0 * h[i] + dt * D1(q, i) - 2.0 * dt * f[i] - 2.0 * dt * ff[i] -
             2.0 * h0[i] + dt * D1(q0, i);

    // Q component
    res[i + N] =
        2.0 * q[i] - 2.0 * dt * f[i + N] - 0.5 * dt * ff[i] * q[i] / h[i] -
        9.0 / 7.0 * dt * q[i] * q[i] / h[i] / h[i] * D1(h, i) +
        5.0 * dt / 3.0 / RE / tan(THETA) * h[i] * D1(h, i) +
        17.0 * dt / 7.0 * q[i] / h[i] * D1(q, i) - 5.0 * dt / 3.0 / RE * h[i] -
        5.0 * dt / 6.0 / CA / RE * h[i] * D3(h, i) +
        5.0 * dt / 2.0 / RE * q[i] / h[i] / h[i] - 2.0 * q0[i] -
        0.5 * dt * ff[i] * q0[i] / h0[i] -
        9.0 / 7.0 * dt * q0[i] * q0[i] / h0[i] / h0[i] * D1(h0, i) +
        5.0 * dt / 3.0 / RE / tan(THETA) * h0[i] * D1(h0, i) +
        17.0 * dt / 7.0 * q0[i] / h0[i] * D1(q0, i) -
        5.0 * dt / 3.0 / RE * h0[i] -
        5.0 * dt / 6.0 / CA / RE * h0[i] * D3(h0, i) +
        5.0 * dt / 2.0 / RE * q0[i] / h0[i] / h0[i];
    res_norm_2 += res[i] * res[i];
    res_norm_2 += res[i + N] * res[i + N];
  }

  return res_norm_2;
}

/* compute the Jacobian matrix for solving the WR implicit time-stepping problem
 */
void est_compute_jacobian(double dt, double **J) {
  double *h = EST_h;
  double *q = EST_q;
  double *ff = EST_ff;

  // TODO: use memset
  for (int i = 0; i < 2 * N; i++) {
    for (int j = 0; j < 2 * N; j++) {
      J[i][j] = 0.0;
    }
  }

  // dFh/dh (top left)
  for (int i = 0; i < N; i++) {
    J[i][i] = 2.0;
  }

  // dFh/dq (top right)
  for (int i = 0; i < N; i++) {
    J[i][N + WRAP(i - 1)] = dt * (-0.5 / DX);
    J[i][N + WRAP(i + 1)] = dt * (0.5 / DX);
  }

  // dFq/dh (bottom left)
  for (int i = 0; i < N; i++) {
    double c1 = dt * (-9.0 / 7.0 * q[i] * q[i] / h[i] / h[i] +
                      5.0 / 3.0 / RE / tan(THETA) * h[i]);
    double c3 = -dt * 5.0 / 6.0 / CA / RE * h[i];
    J[N + i][WRAP(i - 2)] = (-0.5 / DX / DX / DX) * c3;
    J[N + i][WRAP(i - 1)] = (-0.5 / DX) * c1 + (1.0 / DX / DX / DX) * c3;
    J[N + i][WRAP(i + 0)] =
        dt *
        (0.5 * ff[i] * q[i] / h[i] / h[i] +
         18.0 / 7.0 * q[i] * q[i] / h[i] / h[i] / h[i] * D1(h, i) +
         5.0 / 3.0 / RE / tan(THETA) * D1(h, i) -
         17.0 / 7.0 * q[i] / h[i] / h[i] * D1(q, i) - 5.0 / 3.0 / RE -
         5.0 / 6.0 / CA / RE * D3(h, i) - 5.0 / RE * q[i] / h[i] / h[i] / h[i]);
    J[N + i][WRAP(i + 1)] = (0.5 / DX) * c1 + (-1.0 / DX / DX / DX) * c3;
    J[N + i][WRAP(i + 2)] = (0.5 / DX / DX / DX) * c3;
  }

  // dFq/dq (bottom right)
  for (int i = 0; i < N; i++) {
    double c1 = dt * 17.0 / 7.0 * q[i] / h[i];
    J[N + i][N + WRAP(i - 1)] = (-0.5 / DX) * c1;
    J[N + i][N + WRAP(i + 0)] =
        2.0 +
        dt * (-0.5 * ff[i] / h[i] - 18.0 / 7.0 * q[i] / h[i] / h[i] * D1(h, i) +
              17.0 / 7.0 / h[i] * D1(q, i) + 5.0 / 2.0 / RE / h[i] / h[i]);
    J[N + i][N + WRAP(i + 1)] = (0.5 / DX) * c1;
  }
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
    // purely proportional control
    EST_f[i] = 0.0;
    for (int j = 0; j < P; j++) {
      EST_f[i] += EST_L[i][j] * EST_y[j];
    } // j end
  }

  // compute the forcing term for the main system
  for (int i = 0; i < N; i++) {
    EST_ff[i] = control(ITOX(i));
  }

  // implicit time-stepping (for stability)
  // iterate to the solution for the next timestep
  const int iter_max = 100;
  const double res_tol_2 = 1.0e-10; // square of the residual tolerance
  int k = 0;
  for (; k < iter_max; k++) {
    /* compute res */
    double res_norm_2 = est_compute_residual(dt, EST_res);

    /* end if converged */
    if (res_norm_2 < res_tol_2) {
      break;
    }

    /* compute Jacobian and solve linear system */
    // TODO: use sparse/banded matrix representation and solver
    // TODO: since the top block-row of the Jacobian is constant we could
    //       decompose it into a 2x2 block system and save a lot of work
    est_compute_jacobian(dt, EST_J);
    dsv(EST_J, EST_res, 2 * N);

    /* update variables */
    for (int i = 0; i < N; i++) {
      EST_h[i] -= EST_res[i];
    } // i end
    for (int i = 0; i < N; i++) {
      EST_q[i] -= EST_res[N + i];
    } // i end
  } // k end

  /* set estimate at prev time to current time */
  // TODO: pointer swap would be faster
  for (int i = 0; i < N; i++) {
    EST_h0[i] = EST_h[i];
    EST_q0[i] = EST_q[i];
  } // i end

  return k + 1;
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
  EST_J = malloc_f2d(2 * N, 2 * N);

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
  free_2d(EST_J);
}

/* [REQUIRED] steps the system forward in time given the interfacial height */
int est_step(double dt, double *h, int control_on) {
  /* u = K * (h-1) */
  for (int i = 0; i < M; i++) {
    Amag[i] = 0.0;

    if (control_on) {
      for (int j = 0; j < N; j++) {
        Amag[i] += EST_K[i][j] * (EST_h[j] - 1.0);
        Amag[i] += EST_K[i][j + N] * (EST_q[j] - 2.0 / 3.0);
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
