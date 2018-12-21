/*    
    Copyright 2013-2019 Onera.

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

#include "linear.h"

//=============================================================================
/* Produit matrice vecteur
   IN: A: n x m matrix
   IN: B: m vector
   OUT: C: n vector
*/
//=============================================================================
void K_LINEAR::prodv(E_Int n, E_Int m, E_Float* A, E_Float* B, E_Float* C)
{
  E_Float sum = 0.;
  for (E_Int i = 0; i < n; i++)
  {
    sum = 0.;
    for (E_Int j = 0; j < m; j++) sum += A[i+j*n]*B[j];
    C[i] = sum;
  }
}

//=============================================================================
/* Produit matrice matrice
   IN: A: n x m matrix
   IN: B: m x p matrix
   OUT: C: n x p matrix
*/
//=============================================================================
void K_LINEAR::prodm(E_Int n, E_Int m, E_Int p, 
                     E_Float* A, E_Float* B, E_Float* C)
{
  E_Float sum = 0.;
  for (E_Int i = 0; i < n; i++)
  {
    for (E_Int k = 0; k < p; k++)
    {
      sum = 0.;
      for (E_Int j = 0; j < m; j++) sum += A[i+j*n]*B[j+k*m];
      C[i+k*n] = sum;
    }
  }
}
