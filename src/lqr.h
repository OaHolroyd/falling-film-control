#ifndef LQR_H
#define LQR_H


/* ========================================================================== */
/*   FUNCTION DECLARATIONS                                                    */
/* ========================================================================== */
/* computes the optimal gain matrix K for the system (A, B) with quadratic cost
   given by (u, v). n is the system size and m is the control dimension. */
int lqr(double **A, double **B, double u, double v, int n, int m, double **K);


#endif
