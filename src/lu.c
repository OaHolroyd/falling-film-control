#include "lu.h"


/* ========================================================================== */
/*   FUNCTION DEFINITIONS                                                     */
/* ========================================================================== */
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
