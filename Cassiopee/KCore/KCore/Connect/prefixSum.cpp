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
#include "Connect/connect.h"
#include <vector>

// ============================================================================
// Perform an exclusive prefix sum on an array that is a mask comprised solely
// of zeros and ones. Return the total number of ones, that is the total number
// of tagged elements.
// ============================================================================
E_Int K_CONNECT::prefixSum(std::vector<E_Int>& a)
{
  E_Int n = a.size();
  E_Int nthreads = __NUMTHREADS__;
  std::vector<E_Int> threadSums(nthreads+1, 0);

  // Build thread offsets
  #pragma omp parallel
  {
    E_Int ithread = __CURRENT_THREAD__;
    E_Int ai;
    E_Int localSum = 1;  // re-indexing starts at 1

    #pragma omp for schedule(static)
    for (E_Int i = 0; i < n; i++)
    {
      ai = a[i];
      if (ai > 0)
      {
        a[i] = localSum;
        localSum += ai;
      }
    }

    threadSums[ithread+1] = localSum - 1;
  }

  for (E_Int i = 1; i < nthreads+1; i++) threadSums[i] += threadSums[i-1];

  // Add thread offsets to input array
  #pragma omp parallel
  {
    E_Int ithread = __CURRENT_THREAD__;
    E_Int offset = threadSums[ithread];

    #pragma omp for schedule(static)
    for (E_Int i = 0; i < n; i++)
    {
      if (a[i] > 0) a[i] += offset;
    }
  }

  return threadSums[nthreads];  // total number of tagged elements
}
