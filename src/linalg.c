#include "linalg.h"

#include <stdio.h>
#include <math.h>
#include <lapacke.h>

#include "c-utils.h"

#define COMPLEX double complex


/* ========================================================================== */
/*   FUNCTION DEFINITIONS                                                     */
/* ========================================================================== */
/* returns the trace of an n-by-n double matrix A */
double dtr(double **A, int n) {
  double tr = 0.0;

  for (int i = 0; i < n; i++) {
    tr += A[i][i];
  } // i end

  return tr;
}

/* given an n-by-n double matrix A, fills w with the eigenvalues */
void dev(double **A, double complex *w, int n) {
  double *wr = malloc(n*sizeof(double));
  double *wi = malloc(n*sizeof(double));

  /* get real and imaginary parts */
  LAPACKE_dgeev(LAPACK_ROW_MAJOR, 'N', 'N', n, *A, n, wr, wi, NULL, n, NULL, n);

  /* combine */
  for (int i = 0; i < n; i++) {
    w[i] = wr[i] + I*wi[i];
  } // i end

  free(wr);
  free(wi);
}

/* given an n-by-n double matrix A, computes the spectral radius */
void dsr(double **A, int n, double *l) {
  double *wr = malloc(n*sizeof(double));
  double *wi = malloc(n*sizeof(double));

  /* get real and imaginary parts */
  LAPACKE_dgeev(LAPACK_ROW_MAJOR, 'N', 'N', n, *A, n, wr, wi, NULL, n, NULL, n);

  /* find the largest real component */
  *l = wr[0];
  for (int i = 1; i < n; i++) {
    if (*l < wr[i]) { *l = wr[i]; }
  } // i end

  free(wr);
  free(wi);
}

/* solves AX + XA' + Q = 0 (if transpose is 0) or A'X + XA + Q = 0 (otherwise)
   for X, overwriting Q. All matrices are n-by-n */
void dlyap(double **A, double **Q, int n, int transpose) {
  double *wr = malloc(n*sizeof(double));
  double *wi = malloc(n*sizeof(double));

  /* compute Schur form of A */
  int sdim;
  double **Z = malloc_f2d(n, n);
  LAPACKE_dgees(LAPACK_ROW_MAJOR, 'V', 'N', NULL, n, *A, n, &sdim, wr, wi, *Z, n);

  /* compute F = Z'Q */
  double **F = malloc_f2d(n, n);
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      F[i][j] = 0.0;
      for (int k = 0; k < n; k++) {
        F[i][j] += Z[k][i] * Q[k][j];
      } // k end
    } // j end
  } // i end

  // fprintf(stderr, "Q\n  %15lf %15lf\n  %15lf %15lf\n", Q[0][0], Q[0][1], Q[1][0], Q[1][1]);
  // fprintf(stderr, "F\n  %15lf %15lf\n  %15lf %15lf\n", F[0][0], F[0][1], F[1][0], F[1][1]);

  /* compute Q = -FZ */
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      Q[i][j] = 0.0;
      for (int k = 0; k < n; k++) {
        Q[i][j] -= F[i][k] * Z[k][j];
      } // k end
    } // j end
  } // i end

  // fprintf(stderr, "Z\n  %15lf %15lf\n  %15lf %15lf\n", Z[0][0], Z[0][1], Z[1][0], Z[1][1]);
  // fprintf(stderr, "Q\n  %15lf %15lf\n  %15lf %15lf\n", Q[0][0], Q[0][1], Q[1][0], Q[1][1]);
  // ABORT("expected");

  /* solve Schur version */
  double s;
  if (transpose == 0) {
    /* AY + YA' = Q */
    double scale = LAPACKE_dtrsyl(LAPACK_ROW_MAJOR, 'N', 'T', 1, n, n, *A, n, *A, n, *Q, n, &s);
  } else {
    /* A'Y + YA = Q */
    double scale = LAPACKE_dtrsyl(LAPACK_ROW_MAJOR, 'T', 'N', 1, n, n, *A, n, *A, n, *Q, n, &s);
    // fprintf(stderr, "scale: %lf\n", scale);
  }

  s = 1/s;

  /* convert back to original system */
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      F[i][j] = 0.0;
      for (int k = 0; k < n; k++) {
        F[i][j] += Z[i][k] * Q[k][j] * s;
      } // k end
    } // j end
  } // i end
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      Q[i][j] = 0.0;
      for (int k = 0; k < n; k++) {
        Q[i][j] += F[i][k] * Z[j][k];
      } // k end
    } // j end
  } // i end

  free(wr);
  free(wi);
  free_2d(Z);
  free_2d(F);
}


/* ==================== */
/*   LU FACTORISATION   */
/* ==================== */
/* given an n-by-n double matrix A with non-zero determinant, solves Ax = b in
   place. */
void dsv(double **A, double *b, int n) {
  dlu(A, n);
  dlusv(A, b, n);
}

/* given an n-by-n double complex matrix A with non-zero determinant, solves
   Ax = b in place. */
void zsv(double complex **A, double complex *b, int n) {
  zlu(A, n);
  zlusv(A, b, n);
}

/* for an n-by-n double matrix A with non-zero determinant, computes the LU
   factorisation in place with no pivoting */
void dlu(double **A, int n) {
  for (int k = 0; k < n-1; k++) {
    for (int j = k+1; j < n; j++) {
      A[j][k] /= A[k][k];
      for (int i = k+1; i < n; i++) {
        A[j][i] -= A[j][k]*A[k][i];
      } // i end
    } // j end
  } // k end
}

/* for an n-by-n double complex matrix A with non-zero determinant, computes
   the LU factorisation in place with no pivoting */
void zlu(double complex **A, int n) {
  for (int k = 0; k < n-1; k++) {
    for (int j = k+1; j < n; j++) {
      A[j][k] /= A[k][k];
      for (int i = k+1; i < n; i++) {
        A[j][i] -= A[j][k]*A[k][i];
      } // i end
    } // j end
  } // k end
}

/* given an n-by-n double LU factorisation LU with non-zero determinant, solves
   LUx = b in place. */
void dlusv(double **LU, double *b, int n) {
  int i, j;

  /* solve Ly = b */
  double c;
  for (i = 1; i < n; i++) {
    c = 0.0;
    for (j = 1; j < i; j++) {
      c += LU[i][j] * b[j];
    } // j end
    b[i] = (b[i] - c);
  } // i end

  /* solve Uz = y */
  for (i = n-1; i >= 0; i--) {
    c = 0.0;
    for (j = n-1; j >= i+1; j--) {
      c += LU[i][j] * b[j];
    } // j end
    b[i] = (b[i] - c)/LU[i][i];
  } // i end
}

/* given an n-by-n double complex LU factorisation LU with non-zero determinant,
   solves LUz = b in place. */
void zlusv(double complex **LU, double complex *b, int n) {
  int i, j;

  /* solve Ly = b */
  double complex c;
  for (i = 1; i < n; i++) {
    c = 0.0;
    for (j = 1; j < i; j++) {
      c += LU[i][j] * b[j];
    } // j end
    b[i] = (b[i] - c);
  } // i end

  /* solve Uz = y */
  for (i = n-1; i >= 0; i--) {
    c = 0.0;
    for (j = n-1; j >= i+1; j--) {
      c += LU[i][j] * b[j];
    } // j end
    b[i] = (b[i] - c)/LU[i][i];
  } // i end
}


/* =============== */
/*   LQR SOLVERS   */
/* =============== */
/* computes the optimal gain matrix K for the system (A, B) with quadratic cost
   given by (u, v). n is the system size and m is the control dimension. */
   // TODO: add a lqr_work variant for repeated computations
int dlqr(double **A, double **B, double u, double v, int n, int m, double **K) {
  int i, j, k;
  double c = -1/v;

  /* generate the Hamiltonian */
  COMPLEX *H = malloc(2*n*2*n*sizeof(COMPLEX));

  k = 0;
  for (i = 0; i < n; i++) { // top left
    for (j = 0; j < n; j++) {
      H[2*n*i+j] = A[i][j];
    } // j end
  } // i end

  for (i = 0; i < n; i++) { /* top right (-PHI * V^-1 * PSI^T) */
    for (j = i; j < n; j++) {
      H[2*n*i+(j+n)] = 0.0;
      for (k = 0; k < m; k++) {
        H[2*n*i+(j+n)] += B[i][k] * B[j][k];
      } // k end
      H[2*n*i+(j+n)] *= c;
      H[2*n*j+(i+n)] = H[2*n*i+(j+n)]; // use symmetry of PHI * PSI^T
    } // j end
  } // i end

  for (i = 0; i < n; i++) { // bottom left
    for (j = 0; j < n; j++) {
      H[2*n*(i+n)+j] = 0.0;
    } // j end
    H[2*n*(i+n)+i] = -u;
  } // i end

  for (i = 0; i < n; i++) { // bottom right
    for (j = 0; j < n; j++) {
      H[2*n*(i+n)+(j+n)] = -A[j][i];
    } // j end
  } // i end


  /* compute eigenvalues/vectors */
  COMPLEX *w = malloc(2*n*sizeof(COMPLEX));
  COMPLEX *vr = malloc(2*n*2*n*sizeof(COMPLEX));
  int info = LAPACKE_zgeev(LAPACK_ROW_MAJOR, 'N', 'V', 2*n, H,
                           2*n, w, NULL, 2*n, vr, 2*n);


  /* extract stable eigenvectors */ // TODO: remove the need for E
  COMPLEX *E = malloc(2*n*2*n*sizeof(COMPLEX));
  COMPLEX *Q = malloc(n*n*sizeof(COMPLEX));
  COMPLEX *P = malloc(n*n*sizeof(COMPLEX));
  k = 0;
  for (i = 0; i < 2*n; i++) {
    if (creal(w[i]) < 0) {
      for (j = 0; j < 2*n; j++) {
        E[2*n*j+k] = vr[2*n*j+i];
      } // j end
      k++;
    }
  } // i end

  for (i = 0; i < n; i++) {
    for (j = 0; j < n; j++) {
      Q[n*i+j] = E[2*n*i+j];
      P[n*i+j] = E[2*n*(i+n)+j];
    } // j end
  } // i end


  /* compute P */
  int *ipiv = malloc(n*sizeof(int));
  info = LAPACKE_zgesv(LAPACK_COL_MAJOR, n, n, Q, n, ipiv, P, n); // use COL_MAJOR to solve XA=B


  /* compute K */
  for (i = 0; i < m; i++) {
    for (j = 0; j < n; j++) {
      K[i][j] = 0.0;
      for (k = 0; k < n; k++) {
        K[i][j] -= B[k][i] * creal(P[k*n+j]);
      } // k end
      K[i][j] *= c;
    } // j end
  } // i end


  /* free workspace */
  free(H);
  free(E);
  free(Q);
  free(P);
  free(w);
  free(vr);
  free(ipiv);

  return info;
}

/* complex version of the above */
int zlqr(COMPLEX **A, COMPLEX **B, double u, double v, int n, int m, COMPLEX **K) {
  int i, j, k;
  double c = -1/v;

  /* generate the Hamiltonian */
  COMPLEX *H = malloc(2*n*2*n*sizeof(COMPLEX));

  k = 0;
  for (i = 0; i < n; i++) { // top left
    for (j = 0; j < n; j++) {
      H[2*n*i+j] = A[i][j];
    } // j end
  } // i end

  for (i = 0; i < n; i++) { /* top right (-PHI * V^-1 * PSI^T) */
    for (j = i; j < n; j++) {
      H[2*n*i+(j+n)] = 0.0;
      for (k = 0; k < m; k++) {
        H[2*n*i+(j+n)] += B[i][k] * conj(B[j][k]);
      } // k end
      H[2*n*i+(j+n)] *= c;
      H[2*n*j+(i+n)] = H[2*n*i+(j+n)]; // use symmetry of PHI * PSI^T
    } // j end
  } // i end

  for (i = 0; i < n; i++) { // bottom left
    for (j = 0; j < n; j++) {
      H[2*n*(i+n)+j] = 0.0;
    } // j end
    H[2*n*(i+n)+i] = -u;
  } // i end

  for (i = 0; i < n; i++) { // bottom right
    for (j = 0; j < n; j++) {
      H[2*n*(i+n)+(j+n)] = -conj(A[j][i]);
    } // j end
  } // i end


  /* compute eigenvalues/vectors */
  COMPLEX *w = malloc(2*n*sizeof(COMPLEX));
  COMPLEX *vr = malloc(2*n*2*n*sizeof(COMPLEX));
  int info = LAPACKE_zgeev(LAPACK_ROW_MAJOR, 'N', 'V', 2*n, H,
                           2*n, w, NULL, 2*n, vr, 2*n);


  /* extract stable eigenvectors */ // TODO: remove the need for E
  COMPLEX *E = malloc(2*n*2*n*sizeof(COMPLEX));
  COMPLEX *Q = malloc(n*n*sizeof(COMPLEX));
  COMPLEX *P = malloc(n*n*sizeof(COMPLEX));
  k = 0;
  for (i = 0; i < 2*n; i++) {
    if (creal(w[i]) < 0) {
      for (j = 0; j < 2*n; j++) {
        E[2*n*j+k] = vr[2*n*j+i];
      } // j end
      k++;
    }
  } // i end

  for (i = 0; i < n; i++) {
    for (j = 0; j < n; j++) {
      Q[n*i+j] = E[2*n*i+j];
      P[n*i+j] = E[2*n*(i+n)+j];
    } // j end
  } // i end


  /* compute P */
  int *ipiv = malloc(n*sizeof(int));
  info = LAPACKE_zgesv(LAPACK_COL_MAJOR, n, n, Q, n, ipiv, P, n); // use COL_MAJOR to solve XA=B


  /* compute K */
  for (i = 0; i < m; i++) {
    for (j = 0; j < n; j++) {
      K[i][j] = 0.0;
      for (k = 0; k < n; k++) {
        K[i][j] -= conj(B[k][i]) * creal(P[k*n+j]);
      } // k end
      K[i][j] *= c;
    } // j end
  } // i end


  /* free workspace */
  free(H);
  free(E);
  free(Q);
  free(P);
  free(w);
  free(vr);
  free(ipiv);

  return info;
}


/* ========================= */
/*   PENTADIAGONAL SOLVERS   */
/* ========================= */
void pentadiagonal_lu_factorise(const double *a, double *b, double *c, double *d, double *e, int n) {
  // compute the entries of L and U
  d[0] /= c[0];
  e[0] /= c[0];

  c[1] -= b[1] * d[0];
  d[1] = (d[1] - b[1] * e[0]) / c[1];
  e[1] /= c[1];

  for (int i = 2; i < n - 2; i++) {
    b[i] -= a[i] * d[i - 2];
    c[i] -= a[i] * e[i - 2] + b[i] * d[i - 1];
    d[i] = (d[i] - b[i] * e[i - 1]) / c[i];
    e[i] /= c[i];
  }

  b[n - 2] -= a[n - 2] * d[n - 4];
  c[n - 2] -= a[n - 2] * e[n - 4] + b[n - 2] * d[n - 3];
  d[n - 2] = (d[n - 2] - b[n - 2] * e[n - 3]) / c[n - 2];

  b[n - 1] -= a[n - 1] * d[n - 3];
  c[n - 1] -= a[n - 1] * e[n - 3] + b[n - 1] * d[n - 2];
}

void pentadiagonal_lu_solve(const double *al, const double *be, const double *ep, const double *ga, const double *de, double *f, int n) {
  // solve Ly = f via forward substitution
  double *y = f;
  y[0] /= al[0];
  y[1] = (y[1] - be[1] * y[0]) / al[1];
  for (int i = 2; i < n; i++) {
    y[i] = (y[i] - be[i] * y[i - 1] - ep[i] * y[i - 2]) / al[i];
  }

  // solve Ux = y via backward substitution
  double *x = y;
  x[n - 2] -= ga[n - 2] * x[n - 1];
  for (int i = n - 3; i >= 0; i--) {
    x[i] -= ga[i] * x[i + 1] + de[i] * x[i + 2];
  }
}

void pentadiagonal_solve(double *a, double *b, double *c, double *d, double *e, double *f, int n) {
  pentadiagonal_lu_factorise(a, b, c, d, e, n);
  pentadiagonal_lu_solve(c, b, a, d, e, f, n);
}

void cyclic_pentadiagonal_lu_factorise(double *a, double *b, double *c, double *d, double *e, double *k0, double *k1, int n) {
  // set K = [k0 | k1]
  k0[0] = a[0];
  for (int i = 1; i < n - 4; i++) {
    k0[i] = 0.0;
  }
  k0[n - 4] = e[n - 4];
  k0[n - 3] = d[n - 3];

  k1[0] = b[0];
  k1[1] = a[1];
  for (int i = 2; i < n - 3; i++) {
    k1[i] = 0.0;
  }
  k1[n - 3] = e[n - 3];

  // compute the LU factorisation of E
  pentadiagonal_lu_factorise(a, b, c, d, e, n - 2);

  // solve E \ K
  pentadiagonal_lu_solve(c, b, a, d, e, k0, n - 2);
  pentadiagonal_lu_solve(c, b, a, d, e, k1, n - 2);

  // compute the 2x2 matrix C - H E^-1 K for eqn 12
  c[n - 2] -= e[n - 2] * k0[0] + a[n - 2] * k0[n - 4] + b[n - 2] * k0[n - 3];
  d[n - 2] -= e[n - 2] * k1[0] + a[n - 2] * k1[n - 4] + b[n - 2] * k1[n - 3];
  b[n - 1] -= d[n - 1] * k0[0] + e[n - 1] * k0[1] + a[n - 1] * k0[n - 3];
  c[n - 1] -= d[n - 1] * k1[0] + e[n - 1] * k1[1] + a[n - 1] * k1[n - 3];
}

void cyclic_pentadiagonal_lu_solve(double *a, double *b, double *c, double *d, double *e, const double *k0, const double *k1, double *f, int n) {
  // solve E \ f[:-2]
  pentadiagonal_lu_solve(c, b, a, d, e, f, n - 2);

  // compute rhs vector of eqn 12
  f[n - 2] -= e[n - 2] * f[0] + a[n - 2] * f[n - 4] + b[n - 2] * f[n - 3];
  f[n - 1] -= d[n - 1] * f[0] + e[n - 1] * f[1] + a[n - 1] * f[n - 3];

  // solve for the final two elements of the solution vector (eq 12)
  double det = c[n - 2] * c[n - 1] - d[n - 2] * b[n - 1];
  double tmp = (c[n - 1] * f[n - 2] - d[n - 2] * f[n - 1]) / det;
  f[n - 1] = (c[n - 2] * f[n - 1] - b[n - 1] * f[n - 2]) / det;
  f[n - 2] = tmp;

  // update the rhs vector for the final solve (eq 11)
  for (int i = 0; i < n - 2; i++) {
    f[i] -= k0[i] * f[n - 2] + k1[i] * f[n - 1];
  }
}

void cyclic_pentadiagonal_solve(double *a, double *b, double *c, double *d, double *e, double *f, double *k0, double *k1, int n) {
  cyclic_pentadiagonal_lu_factorise(a, b, c, d, e, k0, k1, n);
  cyclic_pentadiagonal_lu_solve(a, b, c, d, e, k0, k1, f, n);
}
