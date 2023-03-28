/*



--------- NUGA v1.0



*/
//Authors : S�m Landier (sam.landier@onera.fr)

#ifndef NUGALIB
#include "parallel.h"

#else

#ifdef _OPENMP

#include <omp.h>
#define __NUMTHREADS__ omp_get_max_threads()
#define __CURRENT_THREAD__  omp_get_thread_num();
#if _OPENMP >= 201307
#define _OPENMP4
#endif

#else // !_OPENMP

#define __NUMTHREADS__ 1
#define __CURRENT_THREAD__ 0

#endif

#define __MIN_SIZE_LIGHT__ 100
#define __MIN_SIZE_MEAN__ 50
#define __MIN_SIZE_HEAVY__ 20

#endif
