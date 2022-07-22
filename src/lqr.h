#ifndef LQR_H
#define LQR_H

#include <complex.h>


/* ========================================================================== */
/*   FUNCTION DECLARATIONS                                                    */
/* ========================================================================== */
/* computes the optimal gain matrix K for the system (A, B) with quadratic cost
   given by (u, v). n is the system size and m is the control dimension. */
int dlqr(double **A, double **B, double u, double v, int n, int m, double **K);

/* complex version of the above */
int zlqr(double complex **A, double complex **B, double u, double v, int n, int m, double complex **K);

int LQR_dgesv(int n, int nrhs, double *a, int lda, int *ipiv, double *b, int ldb);
int LQR_zgesv(int n, int nrhs, double complex *a, int lda, int *ipiv, double complex *b, int ldb);

#endif
