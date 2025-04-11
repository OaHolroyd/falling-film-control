#include "wr.h"

#include <math.h>

#include "linalg.h"

// derivative macros
// NOTE: WRAP requires that 'n' is defined in the current scope, and the
//       derivative macros require that 'dx' is defined in the current scope
#define WRAP(i) ((i + n) % n)

// centre-to-centre spacing
#define D1C(z, i) (0.5 * ((z[WRAP(i + 1)] - z[WRAP(i - 1)]) / dx))

// left and right spacing
#define D0L(z, i) (0.5 * (z[WRAP(i + 1)] + z[WRAP(i)]))
#define D1L(z, i) ((z[WRAP(i + 1)] - z[WRAP(i)]) / dx)
#define D0R(z, i) (0.5 * (z[WRAP(i)] + z[WRAP(i - 1)]))
#define D1R(z, i) ((z[WRAP(i)] - z[WRAP(i - 1)]) / dx)
#define D3R(z, i)                                                              \
  ((-z[WRAP(i - 2)] + 3.0 * z[WRAP(i - 1)] - 3.0 * z[WRAP(i)] +                \
    z[WRAP(i + 1)]) /                                                          \
   (dx * dx * dx))

static double compute_residual(struct wr_data *data, const double dt) {
  // unpack the struct
  const int n = data->n;
  const double dx = data->dx;
  const double theta = data->theta;
  const double re = data->re;
  const double ca = data->ca;
  const double *q = data->q;
  const double *h = data->h;
  const double *h0 = data->h0;
  const double *q0 = data->q0;
  const double *fa = data->fa;
  const double *f = data->f;
  double *resh = data->work;
  double *resq = data->work + n;

  double res_norm_2 = 0.0;

  for (int i = 0; i < n; i++) {
    // cell-centred variables
    const double hc = h[i];
    const double h0c = h0[i];
    const double qxc = D1L(q, i);
    const double q0xc = D1L(q0, i);
    const double fc = f[i];
    const double fac = fa[i];

    // face-centred variables
    const double hf = D0R(h, i);
    const double h0f = D0R(h, i);
    const double hxf = D1R(h, i);
    const double h0xf = D1R(h, i);
    const double hxxxf = D3R(h, i);
    const double h0xxxf = D3R(h, i);
    const double qf = q[i];
    const double q0f = q0[i];
    const double qxf = D1C(q, i);
    const double q0xf = D1C(q0, i);
    const double ff = D0R((f + n), i);
    const double faf = D0R(fa, i);

    // H component (cell-centred)
    resh[i] = hc + 0.5 * dt * qxc - dt * fc - dt * fac - h0c + 0.5 * dt * q0xc;

    // Q component (face-centred)
    resq[i] = qf - dt * ff - 0.25 * dt * faf * qf / hf -
              9.0 / 14.0 * dt * qf * qf / hf / hf * hxf +
              5.0 * dt / 6.0 / re / tan(theta) * hf * hxf +
              17.0 * dt / 14.0 * qf / hf * qxf - 5.0 * dt / 6.0 / re * hf -
              5.0 * dt / 12.0 / ca / re * hf * hxxxf +
              5.0 * dt / 4.0 / re * qf / hf / hf - q0f -
              0.25 * dt * faf * q0f / h0f -
              9.0 / 14.0 * dt * q0f * q0f / h0f / h0f * h0xf +
              5.0 * dt / 6.0 / re / tan(theta) * h0f * h0xf +
              17.0 * dt / 14.0 * q0f / h0f * q0xf - 5.0 * dt / 6.0 / re * h0f -
              5.0 * dt / 12.0 / ca / re * h0f * h0xxxf +
              5.0 * dt / 4.0 / re * q0f / h0f / h0f;
    res_norm_2 += resh[i] * resh[i];
    res_norm_2 += resq[i] * resq[i];
  }

  return res_norm_2;
}

static void newton_step(struct wr_data *data, const double dt) {
  // break J into blocks: J = [[A, B], [C, D]], and note that A = I
  // break res into blocks: res = [a, b]

  // unpack the struct
  const int n = data->n;
  const double dx = data->dx;
  const double theta = data->theta;
  const double re = data->re;
  const double ca = data->ca;
  const double *q = data->q;
  const double *h = data->h;
  const double *fa = data->fa;
  double *a = data->work;
  double *b = data->work + n;
  double *c_l2 = data->work + 2 * n;
  double *c_l1 = data->work + 3 * n;
  double *c_d0 = data->work + 4 * n;
  double *c_u1 = data->work + 5 * n;
  double *s_l2 = data->work + 6 * n;
  double *s_l1 = data->work + 7 * n;
  double *s_d0 = data->work + 8 * n;
  double *s_u1 = data->work + 9 * n;
  double *s_u2 = data->work + 10 * n;
  double *k0 = data->work + 11 * n;
  double *k1 = data->work + 12 * n;
  double *z = data->work + 13 * n;

  const double beta = 1.0 / tan(theta);

  // construct the diagonals of C
  for (int i = 0; i < n; i++) {
    // face-centred variables
    const double hf = D0R(h, i);
    const double hxf = D1R(h, i);
    const double hxxxf = D3R(h, i);
    const double qf = q[i];
    const double qxf = D1C(q, i);
    const double faf = D0R(fa, i);

    // interim constants
    const double c0 = 0.25 * dt * faf * qf / hf / hf +
                      9.0 / 7.0 * dt * qf * qf / hf / hf / hf * hxf +
                      2.5 * dt * beta / 3.0 / re * hxf -
                      17.0 * dt / 14.0 * qf / hf / hf * qxf -
                      5.0 * dt / 6.0 / re - 5.0 * dt / 12.0 / ca / re * hxxxf -
                      2.5 * dt / re * qf / hf / hf / hf;
    const double c1 =
        -9.0 / 14.0 * dt * qf * qf / hf / hf + 5.0 * dt * beta / 6.0 / re * hf;
    const double c3 = -5.0 * dt / 12.0 / ca / re * hf;

    // set diagonals
    c_l2[i] = (-1.0 / dx / dx / dx) * c3;
    c_l1[i] = (3.0 / dx / dx / dx) * c3 + (-1.0 / dx) * c1 + (0.5) * c0;
    c_d0[i] = (-3.0 / dx / dx / dx) * c3 + (1.0 / dx) * c1 + (0.5) * c0;
    c_u1[i] = (1.0 / dx / dx / dx) * c3;
  } // i end

  // construct the diagonals of the Schur complement, D - C A \ B = D - C B
  for (int i = 0; i < n; i++) {
    // face-centred variables
    const double hf = D0R(h, i);
    const double hxf = D1R(h, i);
    const double qf = q[i];
    const double qxf = D1C(q, i);
    const double faf = D0R(fa, i);

    // interim constants
    const double d0 =
        1.0 - 0.25 * dt * faf / hf - 9.0 / 7.0 * dt * qf / hf / hf * hxf +
        17.0 * dt / 14.0 / hf * qxf + 5.0 * dt / 4.0 / re / hf / hf;
    const double d1 = 17.0 * dt / 14.0 * qf / hf;

    // set diagonals
    s_l2[i] = c_l2[i] * 0.5 * dt / dx;
    s_l1[i] = (-0.5 / dx) * d1 + (c_l1[i] - c_l2[i]) * 0.5 * dt / dx;
    s_d0[i] = d0 + (c_d0[i] - c_l1[i]) * 0.5 * dt / dx;
    s_u1[i] = (0.5 / dx) * d1 + (c_u1[i] - c_d0[i]) * 0.5 * dt / dx;
    s_u2[i] = -c_u1[i] * 0.5 * dt / dx;
  } // i end

  // z = C a
  for (int i = 0; i < n; i++) {
    z[i] = a[WRAP(i - 2)] * c_l2[i] + a[WRAP(i - 1)] * c_l1[i] +
           a[WRAP(i + 0)] * c_d0[i] + a[WRAP(i + 1)] * c_u1[i];
  }

  // [b, z] = S \ [b, z]
  cyclic_pentadiagonal_lu_factorise(s_l2, s_l1, s_d0, s_u1, s_u2, k0, k1, n);
  cyclic_pentadiagonal_lu_solve(s_l2, s_l1, s_d0, s_u1, s_u2, k0, k1, b, n);
  cyclic_pentadiagonal_lu_solve(s_l2, s_l1, s_d0, s_u1, s_u2, k0, k1, z, n);

  // b = b - z
  for (int i = 0; i < n; i++) {
    b[i] -= z[i];
  }

  // a = a - B b
  for (int i = 0; i < n; i++) {
    a[i] += (b[i] - b[WRAP(i + 1)]) * 0.5 * dt / dx;
  }
}

int wr_step(struct wr_data *data, const double dt, const double tol,
            const int maxiter) {
  // unpack the struct
  const int n = data->n;
  double *q = data->q;
  double *h = data->h;
  double *h0 = data->h0;
  double *q0 = data->q0;
  double *resh = data->work;
  double *resq = data->work + n;

  int iter = 0;
  for (; iter < maxiter; iter++) {
    // compute residual
    double res_norm_2 = compute_residual(data, dt);

    // finish if converged
    if (sqrt(res_norm_2) < (n * tol) && iter > 0) {
      break;
    }

    // perform Newton iteration step
    newton_step(data, dt);

    // update variables
    for (int i = 0; i < n; i++) {
      h[i] -= resh[i];
    } // i end
    for (int i = 0; i < n; i++) {
      q[i] -= resq[i];
    } // i end
  }

  // store current time as previous time for next step
  for (int i = 0; i < n; i++) {
    h0[i] = h[i];
    q0[i] = q[i];
  } // i end

  return iter + 1;
}
