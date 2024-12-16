#include "c-utils.h"

/* Standard headers */
#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>


/* ========================================================================== */
/*   FUNCTION DEFINITIONS                                                     */
/* ========================================================================== */
/* Allocate memory for a 2D double array and match row indices to the
   corresponding memory locations. */
double** malloc_f2d(int Ni, int Nj) {
  /* allocate row memory */
  double **p_2arr = (double **)malloc(Ni*sizeof(double*));
  if (!p_2arr) { return NULL; }

  /* allocate main memory */
  double *mem = (double *)malloc(Ni*Nj*sizeof(double));
  if (!mem) { free(p_2arr); return NULL; }

  /* match rows to memory */
  for (int i = 0; i < Ni; i++) {
    p_2arr[i] = &(mem[i*Nj]);
  } // i end

  return p_2arr;
}

/* Allocate memory for a 2D complex double array and match row indices to the
   corresponding memory locations. */
COMPLEX** malloc_z2d(int Ni, int Nj) {
  /* allocate row memory */
  // TODO: is this correct
  COMPLEX **p_2arr = (COMPLEX **)malloc(Ni*sizeof(COMPLEX*));
  if (!p_2arr) { return NULL; }

  /* allocate main memory */
  COMPLEX *mem = (COMPLEX *)malloc(Ni*Nj*sizeof(COMPLEX));
  if (!mem) { free(p_2arr); return NULL; }

  /* match rows to memory */
  for (int i = 0; i < Ni; i++) {
    p_2arr[i] = &(mem[i*Nj]);
  } // i end

  return p_2arr;
}

/* Frees memory associated with a 2D array */
void internal_free_2d(void** arr) {
  free(*arr);
  free(arr);
}

/* Aborts with an error message */
void ABORT(const char *format, ...) {
  va_list args;
  va_start(args, format);

  fprintf(stderr, "ERROR: ");
  vfprintf(stderr, format, args);
  fprintf(stderr, "\n");

  va_end(args);

  exit(1);
}

/* output 1D double array to file at fname */
int output_d1d(const char *fname, double *A, int ni) {
  FILE *fp = fopen(fname, "w");

  if (!fp) { return 1; }

  for (int i = 0; i < ni; i++) {
    fprintf(fp, "%lf\n", A[i]);
  } // i end
  fclose(fp);

  return 0;
}

/* output 1D double complex array to files at fname */
int output_z1d(const char *fname, double complex *A, int ni) {
  FILE *fp = fopen(fname, "w");

  if (!fp) { return 1; }

  for (int i = 0; i < ni; i++) {
    fprintf(fp, "%.16lf+%.16lfi\n", creal(A[i]), cimag(A[i]));
  } // i end
  fclose(fp);

  return 0;
}

/* output 2D double array to file at fname */
int output_d2d(const char *fname, double **A, int ni, int nj) {
  FILE *fp = fopen(fname, "w");

  if (!fp) { return 1; }

  for (int i = 0; i < ni; i++) {
    for (int j = 0; j < nj; j++) {
      fprintf(fp, "%.16lf ", A[i][j]);
    } // j end
    fprintf(fp, "\n");
  } // i end
  fclose(fp);

  return 0;
}

/* output 2D double complex array to files at fname */
int output_z2d(const char *fname, double complex **A, int ni, int nj) {
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
