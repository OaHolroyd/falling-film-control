#include "linalg.h"

#include <math.h>
#include <lapacke.h>

#define COMPLEX double complex


/* ========================================================================== */
/*   FUNCTION DEFINITIONS                                                     */
/* ========================================================================== */
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
