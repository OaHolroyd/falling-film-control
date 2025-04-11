#ifndef WR_H
#define WR_H

// wrapper data struct for solving the WR equations
struct wr_data {
  int n; // number of grid points
  double dx; // grid spacing
  double theta; // plate angle
  double re; // Reynolds' number
  double ca; // capillary number

  double *h;   // height, size n
  double *q;   // flux, size n
  double *h0;  // height (previous time step), size n
  double *q0;  // flux (previous time step), size n

  double *fa;  // forcing term (from actuator), size n
  double *f;  // forcing term (arbitrary), size 2 * n

  double *work; // space for time stepping, size 14 * n
};

int wr_step(struct wr_data *data, const double dt, const double tol, const int maxiter);

#endif // WR_H
