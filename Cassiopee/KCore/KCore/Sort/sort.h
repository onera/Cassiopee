/*    
    Copyright 2013-2025 Onera.

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

#ifndef _KCORE_SORT_SORT_H_
#define _KCORE_SORT_SORT_H_

# include "Def/DefTypes.h"
# include "Fld/FldArray.h"
# include "Def/DefFunction.h"

namespace K_SORT
{
/* Sort the array of coordinates (x,y,z) with respect to the norm array, 
   from bound bg to bound bd  
   Quicksort is a recursive algorithm of complexity = O(NlogN) */
  void quickSort(E_Int bg, E_Int bd, K_FLD::FldArrayF& coord, 
                 K_FLD::FldArrayF& norm2);

/* Compute the square norm of any point of coordinates coord(ind,.) */
  void computeNorm(K_FLD::FldArrayF& coord, K_FLD::FldArrayF& norm2);

/* Sort coordinates of points with respect to the square norm using quicksort
   In: coordinates of points
   Out: sorted coordinates with respect to their norm
   corresponding sorted SQUARE norm */
  void sortCoordinates(K_FLD::FldArrayF& coord, K_FLD::FldArrayF& norm2);

/* Make a splitting of the array coord between bg and bd and return 
   the pivot */
  E_Int pivoting(E_Int bg, E_Int bd, K_FLD::FldArrayF& coord, 
                 K_FLD::FldArrayF& norm2);

/* Remove identical points in a sorted vectOfPoints */
  void removeDoublePoints(K_FLD::FldArrayF& vectOfPoints, 
                          K_FLD::FldArrayF& norm2);
}
#endif

// ========================== KCore/Sort/sort.h ========================
