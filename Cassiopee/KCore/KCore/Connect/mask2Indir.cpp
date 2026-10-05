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
// Enumerate the tagged entries of a mask comprised solely of zeros and ones.
// Each entry equal to 1 is replaced by its rank among the tagged entries
// (1-based, in increasing index order); entries equal to 0 are left at 0.
// Example: [0,1,1,0,1] -> [0,1,2,0,3].
// This is an inclusive scan of the mask, evaluated at tagged entries only.
// Return the total number of ones, that is the total number of tagged elements.
// ============================================================================
E_Int K_CONNECT::mask2Indir(std::vector<E_Int>& a)
{
  const E_Int n = a.size();
  const E_Int nthreads = __NUMTHREADS__;
  std::vector<E_Int> threadSums(nthreads+1, 0);

  #pragma omp parallel num_threads(nthreads)
  {
    E_Int ithread = __CURRENT_THREAD__;
    E_Int ai, offset;
    E_Int localSum = 1;  // re-indexing starts at 1

    // Build thread offsets
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

    #pragma omp barrier
    #pragma omp single
    {
      for (E_Int i = 1; i < nthreads+1; i++) threadSums[i] += threadSums[i-1];
    }

    // Add thread offsets to input array
    offset = threadSums[ithread];

    #pragma omp for schedule(static)
    for (E_Int i = 0; i < n; i++)
    {
      if (a[i] > 0) a[i] += offset;
    }
  }

  return threadSums[nthreads];  // total number of tagged elements
}

// ============================================================================
// Enumerate the tagged entries of a mask comprised solely of zeros and ones.
// Each entry equal to 1 is replaced by its rank among the tagged entries of its
// bucket (renumbering restarts at 1 in each bucket).
// This is an inclusive scan of the mask, evaluated at tagged entries only.
// Return the total number of ones per bucket.
// ============================================================================
std::vector<E_Int> K_CONNECT::mask2Indir(std::vector<E_Int>& a,
                                         const std::vector<E_Int>& buckets)
{
  const E_Int n = a.size();
  const E_Int nbuckets = buckets.size();

  std::vector<E_Int> bucketTotals(nbuckets, -1);

  E_Int nn = 0;
  for (E_Int b = 0; b < nbuckets; b++) nn += buckets[b];
  if (nn != n) return bucketTotals;

  const E_Int nthreads = __NUMTHREADS__;
  std::vector<E_Int> threadSums(nbuckets*(nthreads+1));

  // Compute bucket starting offsets
  std::vector<E_Int> bucketStart(nbuckets, 0);
  for (E_Int b = 1; b < nbuckets; b++)
    bucketStart[b] = bucketStart[b-1] + buckets[b-1];

  #pragma omp parallel num_threads(nthreads)
  {
    E_Int ithread = __CURRENT_THREAD__;
    E_Int ai, offset, beg, end;
    E_Int localSum;

    // Build thread offsets per bucket
    for (E_Int b = 0; b < nbuckets; b++)
    {
      localSum = 1;  // re-indexing starts at 1
      beg = bucketStart[b];
      end = beg + buckets[b];

      #pragma omp for nowait schedule(static)
      for (E_Int i = beg; i < end; i++)
      {
        ai = a[i];
        if (ai > 0)
        {
          a[i] = localSum;
          localSum += ai;
        }
      }
      threadSums[(ithread+1)*nbuckets + b] = localSum - 1;
    }

    // Build thread offsets per bucket
    #pragma omp barrier
    #pragma omp single
    {
      for (E_Int b = 0; b < nbuckets; b++)
      {
        for (E_Int i = 0; i < nthreads; i++)
          threadSums[(i+1)*nbuckets + b] += threadSums[i*nbuckets + b];
        bucketTotals[b] = threadSums[nthreads*nbuckets + b];
      }
    }

    // Add thread offsets to input array
    for (E_Int b = 0; b < nbuckets; b++)
    {
      offset = threadSums[ithread*nbuckets + b];
      beg = bucketStart[b];
      end = beg + buckets[b];

      #pragma omp for nowait schedule(static)
      for (E_Int i = beg; i < end; i++)
      {
        if (a[i] > 0) a[i] += offset;
      }
    }
  }

  return bucketTotals;  // total number of tagged elements per buckets
}
