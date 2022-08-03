#ifndef C_UTILS_H
#define C_UTILS_H

#include <complex.h>


#define COMPLEX double complex


/* ========================================================================== */
/*   FUNCTION DECLARATIONS                                                    */
/* ========================================================================== */
/* Allocate memory for a 2D double array and match row indices to the
   corresponding memory locations. */
double** malloc_f2d(int Ni, int Nj);

/* Allocate memory for a 2D complex double array and match row indices to the
   corresponding memory locations. */
COMPLEX** malloc_z2d(int Ni, int Nj);

/* Frees memory associated with a 2D array */
#define free_2d(A) internal_free_2d((void **)A)
void internal_free_2d(void** arr);

/* Aborts with an error message */
void ABORT(const char *format, ...);


#endif
