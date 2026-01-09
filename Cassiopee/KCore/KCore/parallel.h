/*    
    Copyright 2013-2026 ONERA.

    This file is part of Cassiopee.

    Cassiopee is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    Cassiopee is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Cassiopee.  If not, see <http://www.gnu.org/licenses/>.
*/
#ifdef _OPENMP

#include <omp.h>
#define __NUMTHREADS__ omp_get_max_threads()
#define __CURRENT_THREAD__  omp_get_thread_num();

#if _OPENMP >= 201307
#define _OPENMP4
#endif

#else
#define __NUMTHREADS__ 1
#define __CURRENT_THREAD__ 0
#endif

#define __MIN_SIZE_LIGHT__ 100
#define __MIN_SIZE_MEAN__ 50
#define __MIN_SIZE_HEAVY__ 20

#ifndef CACHELINE
#define CACHELINE 32
#endif