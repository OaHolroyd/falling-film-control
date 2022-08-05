#ifndef LINALG_H
#define LINALG_H

#include <complex.h>


/* ========================================================================== */
/*   FUNCTION DECLARATIONS                                                    */
/* ========================================================================== */
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
/* computes the optimal gain matrix K for the system (A, B) with quadratic cost
   given by (u, v). n is the system size and m is the control dimension. */
int dlqr(double **A, double **B, double u, double v, int n, int m, double **K);

/* complex version of the above */
int zlqr(double complex **A, double complex **B, double u, double v, int n, int m, double complex **K);


#endif
