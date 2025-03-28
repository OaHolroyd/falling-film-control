#ifndef CONTROL_ESTIMATOR_H
#define CONTROL_ESTIMATOR_H

#include <math.h>

#include "c-utils.h"
#include "control-core.h"
#include "linalg.h"

#include <string.h>

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

// centre-to-centre spacing
#define D1C(z, i) (0.5 * ((z[WRAP(i + 1)] - z[WRAP(i - 1)]) / DX))

// left and right spacing
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
/* jacobian matrix for WR using face-centred flux */
void wr_jacobian_cf(double **A) {
  memset(A[0], 0, 4 * N * N * sizeof(double));

  // top left (hh): 0

  // top right (hq): -D1L
  const double chq0 = -1.0;
  for (int i = 0; i < N; i++) {
    A[i][N + WRAP(i + 0)] = chq0 * (-1.0 / DX);
    A[i][N + WRAP(i + 1)] = chq0 * (1.0 / DX);
  }

  // bottom left (qh): 5/RE D0R + (4/7-5/3/RE/tan(THETA)) D1R + 5/6/CA/RE D3R
  const double cqh0 = 5.0 / RE;
  const double cqh1 = 4.0 / 7.0 - 5.0 / 3.0 / RE / tan(THETA);
  const double cqh3 = 5.0 / 6.0 / CA / RE;
  for (int i = 0; i < N; i++) {
    A[N + i][WRAP(i - 1)] = cqh3 * (-1.0 / (DX*DX*DX));
    A[N + i][WRAP(i + 0)] = cqh3 * (3.0 / (DX*DX*DX)) + cqh1 * (-1.0 / DX) + cqh0 * (0.5);
    A[N + i][WRAP(i + 1)] = cqh3 * (-3.0 / (DX*DX*DX)) + cqh1 * (1.0 / DX) + cqh0 * (0.5);
    A[N + i][WRAP(i + 2)] = cqh3 * (1.0 / (DX*DX*DX));
  }

  // bottom right (qq): -5/2/RE D0C - 34/21 D1C
  const double cqq0 = -5.0 / 2.0 / RE;
  const double cqq1 = -34.0 / 21.0;
  for (int i = 0; i < N; i++) {
    A[N + i][N + WRAP(i - 1)] = cqq1 * (-0.5 / DX);
    A[N + i][N + WRAP(i + 0)] = cqq0;
    A[N + i][N + WRAP(i + 1)] = cqq1 * (0.5 / DX);
  }
}

/* actuator matrix for WR using face-centred flux */
void wr_actuator_cf(double **B) {
  /* forcing matrix */
  double **F = malloc_f2d(N, M);
  forcing_matrix(F);

  /* actuator matrix */
  for (int i = 0; i < N; i++) {
    for (int j = 0; j < M; j++) {
      B[i][j] = F[i][j];
      // TODO: check this
      B[N+i][j] = (1.0/3.0) * F[i][j];
    } // j end
  } // i end

  free_2d(F);
}

/* solve the (transpose) LQR problem to compute L. This requires EST_C to have
 * been filled with the Benney or WR observer. */
void est_forcing_matrix(double **L) {
  // transpose of L
  double **Lt = malloc_f2d(P, 2 * N);

  // transpose of the system matrix
  double **At = malloc_f2d(2 * N, 2 * N);
  wr_jacobian_cf(At);
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

/* compute the gain matrix mapping the estimator to actuator strengths for the
 * main problem */
void est_gain_matrix(double **K) {
  /* Jacobian */
  double **A = malloc_f2d(2 * N, 2 * N);
  wr_jacobian_cf(A);
  wr_jacobian(A);

  /* actuator matrix */
  double **B = malloc_f2d(2 * N, M);
  wr_actuator_cf(B);

  /* full control matrix */
  dlqr(A, B, DX * MU, 1 - MU, 2 * N, M, K);

  output_d2d("out/Acf.dat", A, 2 * N, 2 * N);
  output_d2d("out/Bcf.dat", B, 2 * N, M);

  free_2d(A);
  free_2d(B);
}

/* compute the residual and return the square of the norm */
double est_compute_residual(double dt, double *res) {
  double res_norm_2 = 0.0;

  for (int i = 0; i < N; i++) {
    // cell-centred variables
    const double hc = EST_h[i];
    const double h0c = EST_h0[i];
    const double qxc = D1L(EST_q, i);
    const double q0xc = D1L(EST_q0, i);
    const double fc = EST_f[i];
    const double ffc = EST_ff[i];

    // face-centred variables
    const double hf = D0R(EST_h, i);
    const double h0f = D0R(EST_h, i);
    const double hxf = D1R(EST_h, i);
    const double h0xf = D1R(EST_h, i);
    const double hxxxf = D3R(EST_h, i);
    const double h0xxxf = D3R(EST_h, i);
    const double qf = EST_q[i];
    const double q0f = EST_q0[i];
    const double qxf = D1C(EST_q, i);
    const double q0xf = D1C(EST_q0, i);
    const double ff = EST_f[i + N];
    const double fff = D0R(EST_ff, i);

    // H component (cell-centred)
    res[i] = 2.0 * hc + dt * qxc - 2.0 * dt * fc - 2.0 * dt * ffc - 2.0 * h0c +
             dt * q0xc;

    // Q component (face-centred)
    res[i + N] = 2.0 * qf - 2.0 * dt * ff - 0.5 * dt * fff * qf / hf -
                 9.0 / 7.0 * dt * qf * qf / hf / hf * hxf +
                 5.0 * dt / 3.0 / RE / tan(THETA) * hf * hxf +
                 17.0 * dt / 7.0 * qf / hf * qxf - 5.0 * dt / 3.0 / RE * hf -
                 5.0 * dt / 6.0 / CA / RE * hf * hxxxf +
                 5.0 * dt / 2.0 / RE * qf / hf / hf - 2.0 * q0f -
                 0.5 * dt * fff * q0f / h0f -
                 9.0 / 7.0 * dt * q0f * q0f / h0f / h0f * h0xf +
                 5.0 * dt / 3.0 / RE / tan(THETA) * h0f * h0xf +
                 17.0 * dt / 7.0 * q0f / h0f * q0xf -
                 5.0 * dt / 3.0 / RE * h0f -
                 5.0 * dt / 6.0 / CA / RE * h0f * h0xxxf +
                 5.0 * dt / 2.0 / RE * q0f / h0f / h0f;
    res_norm_2 += res[i] * res[i];
    res_norm_2 += res[i + N] * res[i + N];
  }

  return res_norm_2;
}

/* compute the Jacobian matrix for solving the WR implicit time-stepping problem
 */
void est_compute_jacobian(double dt, double **J) {
  const double BETA = 1.0 / tan(THETA);

  // TODO: use memset
  memset(J[0], 0, 4 * N * N * sizeof(double));

  // dFh/dh (top left)
  for (int i = 0; i < N; i++) {
    J[i][i] = 2.0;
  }

  // dFh/dq (top right)
  for (int i = 0; i < N; i++) {
    const double c0 = dt;
    J[i][N + WRAP(i + 0)] = (-1.0 / DX) * c0;
    J[i][N + WRAP(i + 1)] = (1.0 / DX) * c0;
  }

  // dFq/dh (bottom left)
  for (int i = 0; i < N; i++) {
    // face-centred variables
    const double hf = D0R(EST_h, i);
    const double hxf = D1R(EST_h, i);
    const double hxxxf = D3R(EST_h, i);
    const double qf = EST_q[i];
    const double qxf = D1C(EST_q, i);
    const double fff = D0R(EST_ff, i);

    const double c0 = 0.5 * dt * fff * qf / hf / hf +
                      18.0 / 7.0 * dt * qf * qf / hf / hf / hf * hxf +
                      5.0 * dt * BETA / 3.0 / RE * hxf -
                      17.0 * dt / 7.0 * qf / hf / hf * qxf -
                      5.0 * dt / 3.0 / RE - 5.0 * dt / 6.0 / CA / RE * hxxxf -
                      5.0 * dt / RE * qf / hf / hf / hf;
    const double c1 =
        -9.0 / 7.0 * dt * qf * qf / hf / hf + 5.0 * dt * BETA / 3.0 / RE * hf;
    const double c3 = -5.0 * dt / 6.0 / CA / RE * hf;

    J[N + i][WRAP(i - 1)] += (0.5) * c0;
    J[N + i][WRAP(i + 0)] += (0.5) * c0;

    J[N + i][WRAP(i - 1)] += (-1.0 / DX) * c1;
    J[N + i][WRAP(i + 0)] += (1.0 / DX) * c1;

    J[N + i][WRAP(i - 2)] = (-1.0 / DX / DX / DX) * c3;
    J[N + i][WRAP(i - 1)] = (3.0 / DX / DX / DX) * c3;
    J[N + i][WRAP(i + 0)] = (-3.0 / DX / DX / DX) * c3;
    J[N + i][WRAP(i + 1)] = (1.0 / DX / DX / DX) * c3;
  }

  // dFq/dq (bottom right)
  for (int i = 0; i < N; i++) {
    // face-centred variables
    const double hf = D0R(EST_h, i);
    const double hxf = D1R(EST_h, i);
    const double qf = EST_q[i];
    const double qxf = D1C(EST_q, i);
    const double fff = D0R(EST_ff, i);

    const double c0 =
        2.0 - 0.5 * dt * fff / hf - 18.0 / 7.0 * dt * qf / hf / hf * hxf +
        17.0 * dt / 7.0 / hf * qxf + 5.0 * dt / 2.0 / RE / hf / hf;
    const double c1 = 17.0 * dt / 7.0 * qf / hf;

    J[N + i][N + WRAP(i - 1)] = (-0.5 / DX) * c1;
    J[N + i][N + WRAP(i + 0)] = c0;
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
  const double res_tol = 1.0e-6; // square of the residual tolerance
  int k = 0;
  for (; k < iter_max; k++) {
    /* compute res */
    double res_norm_2 = est_compute_residual(dt, EST_res);

    /* end if converged */
    if (sqrt(res_norm_2) < (N * res_tol) && k > 0) {
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
