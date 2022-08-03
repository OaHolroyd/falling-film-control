#ifndef LU_H
#define LU_H

#include <complex.h>


/* ========================================================================== */
/*   FUNCTION DECLARATIONS                                                    */
/* ========================================================================== */
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


#endif
