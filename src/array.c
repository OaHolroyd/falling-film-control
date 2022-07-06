#include "array.h"

/* Standard headers */
#include <stdlib.h>


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


/* Frees memory associated with a 2D array */
void free_2d(void** arr) {
  free(*arr);
  free(arr);
}
