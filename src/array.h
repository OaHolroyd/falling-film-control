#ifndef ARRAY_H
#define ARRAY_H


/* ========================================================================== */
/*   FUNCTION DECLARATIONS                                                    */
/* ========================================================================== */
/* Allocate memory for a 2D double array and match row indices to the
   corresponding memory locations. */
double** malloc_f2d(int Ni, int Nj);


/* Frees memory associated with a 2D array */
void free_2d(void** arr);


#endif
