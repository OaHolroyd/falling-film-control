#ifndef PARALLEL_H
#define PARALLEL_H

/* openmp seems to be slower on macOS so don't use unless flagged */
#ifdef PARALLEL
#include <omp.h>
#endif

#endif
