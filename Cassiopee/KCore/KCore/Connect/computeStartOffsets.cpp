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

// ============================================================================
// Given a sorted array of global indices belonging to consecutive buckets,
// compute the starting offset in the array for each bucket using binary search.
// IN: sortedIndices: sorted global indices
// IN: nindices is the total number of global indices
// IN: bucketSizes are the number of global indices present in each bucket
// OUT: offsets is of size nbuckets+1, with offsets[0] = 0, and where
// offsets[b] to offsets[b+1]-1 is the range of items in bucket b.
// If a bucket has no items, offsets[b+1] == offsets[b].
// ============================================================================
void K_CONNECT::computeStartOffsets(
  const E_Int* sortedIndices, E_Int nindices,
  const std::vector<E_Int>& bucketSizes,
  std::vector<E_Int>& offsets
)
{
  E_Int lo, mid, hi;
  E_Int nbuckets = bucketSizes.size();
  offsets.clear(); offsets.resize(nbuckets+1); 
  offsets[0] = 0;

  E_Int boundary = 1;
  for (E_Int b = 0; b < nbuckets; b++)
  {
    boundary += bucketSizes[b];
    lo = offsets[b];
    hi = nindices;
    while (lo < hi)
    {
      mid = (lo + hi)/2;
      if (sortedIndices[mid] < boundary) lo = mid + 1;
      else hi = mid;
    }
    offsets[b+1] = lo;
  }
}

