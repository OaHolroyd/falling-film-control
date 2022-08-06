#ifndef DEBUG_H
#define DEBUG_H

#include <stdio.h>
#include <complex.h>


/* ========================================================================== */
/*   FUNCTION DEFINITIONS                                                     */
/* ========================================================================== */
/* output 1D double array to file at fname */
int debug_out_d1d(const char *fname, double *A, int ni) {
  FILE *fp = fopen(fname, "w");

  if (!fp) { return 1; }

  for (int i = 0; i < ni; i++) {
    fprintf(fp, "%lf\n", A[i]);
  } // i end
  fclose(fp);

  return 0;
}

/* output 1D double complex array to files at fname */
int debug_out_z1d(const char *fname, double complex *A, int ni) {
  FILE *fp = fopen(fname, "w");

  if (!fp) { return 1; }

  for (int i = 0; i < ni; i++) {
    fprintf(fp, "%lf+%lfi\n", creal(A[i]), cimag(A[i]));
  } // i end
  fclose(fp);

  return 0;
}

/* output 2D double array to file at fname */
int debug_out_d2d(const char *fname, double **A, int ni, int nj) {
  FILE *fp = fopen(fname, "w");

  if (!fp) { return 1; }

  for (int i = 0; i < ni; i++) {
    for (int j = 0; j < nj; j++) {
      fprintf(fp, "%lf ", A[i][j]);
    } // j end
    fprintf(fp, "\n");
  } // i end
  fclose(fp);

  return 0;
}

/* output 2D double complex array to files at fname */
int debug_out_z2d(const char *fname, double complex **A, int ni, int nj) {
  FILE *fp = fopen(fname, "w");

  if (!fp) { return 1; }

  for (int i = 0; i < ni; i++) {
    for (int j = 0; j < nj; j++) {
      fprintf(fp, "%lf+%lfi ", creal(A[i][j]), cimag(A[i][j]));
    } // j end
    fprintf(fp, "\n");
  } // i end
  fclose(fp);

  return 0;
}


#endif
