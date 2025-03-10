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

# include <stdio.h>
# include "Sort/sort.h"
using namespace K_FLD; 

//=============================================================================
/* Sort the array of coordinates (x,y,z) with respect to the norm array, 
   from bound bg to bound bd  
   Quicksort is a recursive algorithm */
//=============================================================================
void K_SORT::quickSort(E_Int bg, E_Int bd, FldArrayF& coord, FldArrayF& norm2)
{
  if (bg < bd)
  { 
    E_Int ipiv = pivoting(bg, bd, coord, norm2);
    quickSort(bg, ipiv-1, coord, norm2);
    quickSort(ipiv+1, bd, coord, norm2);
  }
}
//=============================================================================
/* 
   Compute the square norm of any point of coordinates coord(ind,.)
   Cette routine malloc norm2.
 */
//=============================================================================
void K_SORT::computeNorm(FldArrayF& coord, FldArrayF& norm2)
{
  E_Int nPts = coord.getSize();
  norm2.malloc(nPts);
  norm2.setAllValuesAtNull();
  E_Float* norm2p = norm2.begin();

  for (E_Int dim = 1; dim <= 3; dim++)
  {
    E_Float* coordp = coord.begin(dim);
    for (E_Int ind = 0; ind < nPts; ind++)
      norm2p[ind] += coordp[ind] * coordp[ind];
  }
}

//=============================================================================
/* Sort coordinates of points with respect to the square norm using quicksort*/
//=============================================================================
void K_SORT::sortCoordinates(FldArrayF& coord, FldArrayF& norm2)
{
  /* Compute the square norm of coord */
  computeNorm(coord, norm2);
 
  /* Sort the coord array with respect to the square norm */
  E_Int bg = 0;
  E_Int bd = coord.getSize()-1;
  quickSort(bg, bd, coord, norm2);
}

//=============================================================================
/* Make a splitting of the array coord between bg and bd and return the pivot*/
//=============================================================================
E_Int K_SORT::pivoting(E_Int bg, E_Int bd, FldArrayF& coord, FldArrayF& norm2)
{
  E_Int nfld = coord.getNfld();
  E_Float* norm2p = norm2.begin();
  E_Float piv = norm2p[bg];
  E_Int l = bg+1;
  E_Int r = bd;
  
  E_Float tmp;

  while (l <= r)
  {
    while (l <= bd && norm2p[l] <= piv) l++;
    while (r >=  0 && norm2p[r] > piv) r--;
    if ( l < r )
    {
      // swap norm(indr) and norm(indl)
      tmp = norm2p[r];
      norm2p[r] = norm2p[l];
      norm2p[l] = tmp;

      // swap indl and indr elements of coord
      for (E_Int eq = 1; eq <= nfld; eq++)
      {
        E_Float* coord1 = coord.begin(eq);
        tmp = coord1[r];
        coord1[r] = coord1[l];
        coord1[l] = tmp;
      }
      r--;
      l++;
    }
  }
  
  // swap norm(indr) and norm(bg)
  tmp = norm2p[r];
  norm2p[r] = norm2p[bg];
  norm2p[bg] = tmp;
  
  // swap indl and bg elements of coord
  for (E_Int eq = 1; eq <= nfld; eq++)
  {
    E_Float* coord1 = coord.begin(eq);
    tmp = coord1[r];
    coord1[r] = coord1[bg];
    coord1[bg] = tmp;
  }
  return r;
}

//=============================================================================
/* Remove identical points in a sorted vectOfPoints */
//=============================================================================
void K_SORT::removeDoublePoints(FldArrayF& vectOfPoints, FldArrayF& norm2)
{
  E_Int npts = vectOfPoints.getSize();
  E_Int nfld = vectOfPoints.getNfld();
  FldArrayF work(npts, nfld);
  FldArrayF workn2(npts);
  E_Float dx, dy, dz;
  E_Boolean found;
  E_Int c=0;
  E_Float* vx = vectOfPoints.begin(1);
  E_Float* vy = vectOfPoints.begin(2);
  E_Float* vz = vectOfPoints.begin(3);
  E_Float* norm2p = norm2.begin();
  E_Float* workn2p = workn2.begin();
  E_Int j;

  for (E_Int i = 0; i < npts; i++)
  {
    found = false;
   
    j = i-1; 
    while (j >=0 && 
           K_FUNC::fEqualZero(norm2p[i]-norm2p[j], 1.e-6) == true)
    {
      dx = vx[i] - vx[j];
      dy = vy[i] - vy[j];
      dz = vz[i] - vz[j];

      if (dx*dx + dy*dy + dz*dz < 1.e-12)
      {
        found = true; break;
      }
      j--;
    }
    if (found == false)
    {
      for (E_Int eq = 1; eq <= nfld; eq++)
      {
        work(c, eq) = vectOfPoints(i, eq);
      }
      workn2p[c] = norm2p[i];
      c++;
    }
  }
  
  work.reAllocMat(c, nfld);
  workn2.reAlloc(c);
  vectOfPoints = work;
  norm2 = workn2;
}
//---------------------- KCore/Sort/sort.cpp -------------------------------
