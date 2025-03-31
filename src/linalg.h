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
/* computes the optimal gain matrix K for the system (A, B) with quadratic cost
   given by (u, v). n is the system size and m is the control dimension. */
int dlqr(double **A, double **B, double u, double v, int n, int m, double **K);

/* complex version of the above */
int zlqr(double complex **A, double complex **B, double u, double v, int n, int m, double complex **K);


/* ========================= */
/*   PENTADIAGONAL SOLVERS   */
/* ========================= */
/**
 * Factorises a pentadiagonal matrix A = pdiag(a, b, c, d, e) into A = LU
 *
 * @param a second lower diagonal band (OVERWRITTEN with the main diagonal band of L)
 * @param b first lower diagonal band (OVERWRITTEN with the first lower diagonal band of L)
 * @param c main diagonal band (OVERWRITTEN with the second lower diagonal band of L)
 * @param d first upper diagonal band (OVERWRITTEN with first upper diagonal band of U)
 * @param e second upper diagonal band (OVERWRITTEN with second upper diagonal band of U)
 * @param n size of the matrix (square)
 */
void pentadiagonal_lu_factorise(const double *a, double *b, double *c, double *d, double *e, int n);

/**
 * Given an LU factorised pentadiagonal matrix A = LU, solves Ax = f in place.
 *
 * @param al main diagonal band of L
 * @param be first lower diagonal band of L
 * @param ep second lower diagonal band of L
 * @param ga first upper diagonal band of U
 * @param de second upper diagonal band of U
 * @param f right-hand side vector (OVERWRITTEN with the solution)
 * @param n size of the matrix (square)
 */
void pentadiagonal_lu_solve(const double *al, const double *be, const double *ep, const double *ga, const double *de, double *f, int n);

/**
 * Solves the system Ax = f in place, where A is a pentadiagonal.
 *
 * @param a second lower diagonal band (OVERWRITTEN with the main diagonal band of L)
 * @param b first lower diagonal band (OVERWRITTEN with the first lower diagonal band of L)
 * @param c main diagonal band (OVERWRITTEN with the second lower diagonal band of L)
 * @param d first upper diagonal band (OVERWRITTEN with first upper diagonal band of U)
 * @param e second upper diagonal band (OVERWRITTEN with second upper diagonal band of U)
 * @param f right-hand side vector (OVERWRITTEN with the solution)
 * @param n size of the matrix (square)
 */
void pentadiagonal_solve(double *a, double *b, double *c, double *d, double *e, double *f, int n);

/**
 * Compute the subsystem LU factorisation and prepares a periodic pentadiagonal matrix for solving.
 *
 * See `cyclic_pentadiagonal_solve` for notation.
 *
 * @param a second lower diagonal band (OVERWRITTEN with the main diagonal band of E)
 * @param b first lower diagonal band (OVERWRITTEN with the firs lower diagonal band of E)
 * @param c main diagonal band (OVERWRITTEN with the second lower diagonal band of E)
 * @param d first upper diagonal band (OVERWRITTEN with the first upper diagonal band of E)
 * @param e second upper diagonal band (OVERWRITTEN with the second upper diagonal band of E)
 * @param k0 OVERWRITTEN with the first column of E^-1 K
 * @param k1 OVERWRITTEN with the second column of E^-1 K
 * @param n size of the matrix (square)
 */
void cyclic_pentadiagonal_lu_factorise(double *a, double *b, double *c, double *d, double *e, double *k0, double *k1, int n);

/**
 * Starting with a system proccesed by `cyclic_pentadiagonal_lu_factorise`,
 * solve Ax=f.
 *
 * @param a the main diagonal band of E
 * @param b the first lower diagonal band of E
 * @param c the second lower diagonal band of E
 * @param d the first upper diagonal band of E
 * @param e the second upper diagonal band of E
 * @param k0 the first column of E^-1 K
 * @param k1 the second column of E^-1 K
 * @param f the right-hand side vector (OVERWRITTEN with the solution)
 * @param n size of the matrix (square)
 */
void cyclic_pentadiagonal_lu_solve(double *a, double *b, double *c, double *d, double *e, const double *k0, const double *k1, double *f, int n);

/**
 * Solve a periodic pentadiagonal system of equations.
 *
 * @param a second lower diagonal band (OVERWRITTEN with the main diagonal band
 * of E)
 * @param b first lower diagonal band (OVERWRITTEN with the firs lower diagonal
 * band of E)
 * @param c main diagonal band (OVERWRITTEN with the second lower diagonal band
 * of E)
 * @param d first upper diagonal band (OVERWRITTEN with the first upper diagonal
 * band of E)
 * @param e second upper diagonal band (OVERWRITTEN with the second upper
 * diagonal band of E)
 * @param f the right-hand side vector (OVERWRITTEN with the solution)
 * @param k0 working space for the first column of E^-1 K
 * @param k1 working space for the second column of E^-1 K
 * @param n size of the matrix (square)
 */
void cyclic_pentadiagonal_solve(double *a, double *b, double *c, double *d, double *e, double *f, double *k0, double *k1, int n);


#endif
