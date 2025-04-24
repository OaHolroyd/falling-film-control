#ifndef BENNEY_H
#define BENNEY_H

// wrapper data struct for solving the Benney equation
struct benney_data {
  int n; // number of grid points
  double dx; // grid spacing
  double theta; // plate angle
  double re; // Reynolds' number
  double ca; // capillary number

  double *h;   // height, size n
  double *h0;  // height (previous time step), size n

  double *fa;  // forcing term (from actuator), size n
  double *f;  // forcing term (arbitrary), size n

  double *work; // space for time stepping, size 8 * n
};

int benney_step(struct benney_data *data, const double dt, const double tol, const int maxiter);

#endif // BENNEY_H
