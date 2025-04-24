#include "benney.h"

#include <math.h>

#include "linalg.h"

// derivative macros
// NOTE: WRAP requires that 'n' is defined in the current scope, and the
//       derivative macros require that 'dx' is defined in the current scope
#define WRAP(i) ((i + n) % n)

// centre-to-centre spacing
#define D1(z, i) (0.5 * ((z[WRAP(i + 1)] - z[WRAP(i - 1)]) / dx))
#define D2(z, i)                                                               \
  (((z[WRAP(i + 1)] - 2.0 * z[WRAP(i)] + z[WRAP(i - 1)]) / (dx * dx)))
#define D3(z, i)                                                               \
  (0.5 *                                                                       \
   (z[WRAP(i + 2)] - 2.0 * z[WRAP(i + 1)] + 2.0 * z[WRAP(i - 1)] -             \
    z[WRAP(i - 2)]) /                                                          \
   (dx * dx * dx))
#define D4(z, i)                                                               \
  ((z[WRAP(i + 2)] - 4.0 * z[WRAP(i + 1)] + 6.0 * z[WRAP(i)] -                 \
    4.0 * z[WRAP(i - 1)] + z[WRAP(i - 2)]) /                                   \
   (dx * dx * dx * dx))

static double compute_residual(struct benney_data *data, const double dt) {
  // unpack the struct
  const int n = data->n;
  const double dx = data->dx;
  const double theta = data->theta;
  const double re = data->re;
  const double ca = data->ca;
  double *res = data->work;

  double res_norm_2 = 0.0;

  for (int i = 0; i < n; i++) {
    const double h = data->h[i];
    const double hx = D1(data->h, i);
    const double hxx = D2(data->h, i);
    const double hxxx = D3(data->h, i);
    const double hxxxx = D4(data->h, i);
    const double h0 = data->h0[i];
    const double h0x = D1(data->h0, i);
    const double h0xx = D2(data->h0, i);
    const double h0xxx = D3(data->h0, i);
    const double h0xxxx = D4(data->h0, i);
    const double fa = data->fa[i];
    const double fax = D1(data->fa, i);
    const double f = data->f[i];

    // H component
    const double H =
        f + fa - 2.0 * h * h * hx + 2.0 / tan(theta) * h * h * hx * hx +
        2.0 / 3.0 / tan(theta) * h * h * h * hxx -
        2.0 / ca * h * h * hx * hxxx - 2.0 / 3.0 / ca * h * h * h * hxxxx -
        16.0 * re / 5.0 * h * h * h * h * h * hx * hx -
        8.0 * re / 15.0 * h * h * h * h * h * h * hxx +
        8.0 * re / 3.0 * h * h * h * hx * fa +
        2.0 * re / 3.0 * h * h * h * h * fax;
    const double H0 = f + fa - 2.0 * h0 * h0 * h0x +
                      2.0 / tan(theta) * h0 * h0 * h0x * h0x +
                      2.0 / 3.0 / tan(theta) * h0 * h0 * h0 * h0xx -
                      2.0 / ca * h0 * h0 * h0x * h0xxx -
                      2.0 / 3.0 / ca * h0 * h0 * h0 * h0xxxx -
                      16.0 * re / 5.0 * h0 * h0 * h0 * h0 * h0 * h0x * h0x -
                      8.0 * re / 15.0 * h0 * h0 * h0 * h0 * h0 * h0 * h0xx +
                      8.0 * re / 3.0 * h0 * h0 * h0 * h0x * fa +
                      2.0 * re / 3.0 * h0 * h0 * h0 * h0 * fax;
    res[i] = 2.0 * h - dt * H - 2.0 * h0 - dt * H0;

    res_norm_2 += res[i] * res[i];
  }

  return res_norm_2;
}

static void newton_step(struct benney_data *data, const double dt) {
  // unpack the struct
  const int n = data->n;
  const double dx = data->dx;
  const double theta = data->theta;
  const double re = data->re;
  const double ca = data->ca;
  double *res = data->work;
  double *j_l2 = data->work + n;
  double *j_l1 = data->work + 2 * n;
  double *j_d0 = data->work + 3 * n;
  double *j_u1 = data->work + 4 * n;
  double *j_u2= data->work + 5 * n;
  double *k0 = data->work + 6 * n;
  double *k1 = data->work + 7 * n;

  const double beta = 1.0 / tan(theta);

  // construct the diagonals of J
  for (int i = 0; i < n; i++) {
    const double h = data->h[i];
    const double hx = D1(data->h, i);
    const double hxx = D2(data->h, i);
    const double hxxx = D3(data->h, i);
    const double hxxxx = D4(data->h, i);
    const double fa = data->fa[i];
    const double fax = D1(data->fa, i);

    // interim constants
    const double c0 = 2.0 - dt * (-4.0*h*hx+4.0*beta*h*hx*hx+2.0*beta*h*h*hxx-4.0/ca*h*hx*hxxx-2.0/ca*h*h*hxxxx-16.0*re*h*h*h*h*hx*hx-16.0*re/5.0*h*h*h*h*h*hxx+8.0*re*h*h*hx*fa + 8.0*re/3.0*h*h*h*fax);
    const double c1 = -dt*(2.0*h*h+4.0*beta*h*h*hx-2.0/ca*h*h*hxxx-32.0*re/5.0*h*h*h*h*h*hx + 8.0*re/3.0*h*h*h*fa);
    const double c2 = -dt * (2.0*beta/3.0*h*h*h - 8.0*re/15.0*h*h*h*h*h*h);
    const double c3 = -dt * (-2.0/ca*h*h*hx);
    const double c4 = -dt * (-2.0/3.0/ca*h*h*h);

    // set diagonals
    j_l2[i] = 1.0/(dx*dx*dx*dx)*c4 - 0.5/(dx*dx*dx)*c3;
    j_l1[i] = -4.0/(dx*dx*dx*dx)*c4 + 1.0/(dx*dx*dx)*c3 + 1.0/(dx*dx)*c2 - 0.5/dx*c1;
    j_d0[i] = 6.0/(dx*dx*dx*dx)*c4 - 2.0/(dx*dx)*c2 + c0;
    j_u1[i] = -4.0/(dx*dx*dx*dx)*c4 - 1.0/(dx*dx*dx)*c3 + 1.0/(dx*dx)*c2 + 0.5/dx*c1;
    j_u2[i] = 1.0/(dx*dx*dx*dx)*c4 + 0.4/(dx*dx*dx)*c3;
  } // i end

  cyclic_pentadiagonal_solve(j_l2, j_l1, j_d0, j_u1, j_u2, res, k0, k1, n);
}

int benney_step(struct benney_data *data, const double dt, const double tol,
                const int maxiter) {
  // unpack the struct
  const int n = data->n;
  double *h = data->h;
  double *h0 = data->h0;
  double *res = data->work;

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
      h[i] -= res[i];
    } // i end
  }

  // store current time as previous time for next step
  for (int i = 0; i < n; i++) {
    h0[i] = h[i];
  } // i end

  return iter + 1;
}
