#ifndef LINALG_H
#define LINALG_H

#include <complex.h>


/* ========================================================================== */
/*   FUNCTION DECLARATIONS                                                    */
/* ========================================================================== */
/* returns the trace of an n-by-n double matrix A */
double dtr(double **A, int n);

/* given an n-by-n double matrix A, fills w with the eigenvalues */
void dev(double **A, double complex *w, int n);

/* given an n-by-n double matrix A, computes the spectral radius l */
void dsr(double **A, int n, double *l);

/* solves AX + XA' + Q = 0 (if transpose is 0) or A'X + XA + Q = 0 (otherwise)
   for X, overwriting Q. All matrices are n-by-n */
void dlyap(double **A, double **Q, int n, int transpose);


/* ==================== */
/*   LU FACTORISATION   */
/* ==================== */
/* given an n-by-n double matrix A with non-zero determinant, solves Ax = b in
   place. */
void dsv(double **A, double *b, int n);

/* given an n-by-n double complex matrix A with non-zero determinant, solves
   Ax = b in place. */
void zsv(double complex **A, double complex *b, int n);

/* for an n-by-n double matrix A with non-zero determinant, computes the LU
   factorisation in place with no pivoting */
void dlu(double **A, int n);

/* for an n-by-n double complex matrix A with non-zero determinant, computes
   the LU factorisation in place with no pivoting */
void zlu(double complex **A, int n);

/* given an n-by-n double LU factorisation LU with non-zero determinant, solves
   LUx = b in place. */
void dlusv(double **LU, double *b, int n);

/* given an n-by-n double complex LU factorisation LU with non-zero determinant,
   solves LUz = b in place. */
void zlusv(double complex **LU, double complex *b, int n);


/* =============== */
/*   LQR SOLVERS   */
/* =============== */
/* solves the Riccati problem associated with (A, B, u, v), filling P with the
   solution. P have space for n*n COMPLEXes. */
int riccati(double **A, double **B, double u, double v, int n, int m, double complex *P);

/* computes the optimal gain matrix K for the system (A, B) with quadratic cost
   given by (u, v). n is the system size and m is the control dimension. */
int dlqr(double **A, double **B, double u, double v, int n, int m, double **K);

/* complex version of the above */
int zlqr(double complex **A, double complex **B, double u, double v, int n, int m, double complex **K);


#endif
