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

#include "linear.h"

//=============================================================================
/* Inverse a 2x2 matrix.
   IN: A: 2x2 matrix
   OUT: B: 2x2 matrix
   Retourne 0 si la matrice n'est pas inversible */
//=============================================================================
E_Int K_LINEAR::inv2(E_Float* A, E_Float* B)
{
  E_Float det;
  E_Float a11 = A[0]; E_Float a21 = A[1]; 
  E_Float a12 = A[2]; E_Float a22 = A[3];

  det = a11*a22 - a12*a21;
  if (K_FUNC::E_abs(det) < 1.e-12) return 0;

  det = 1./det;
  B[0] = det*a22;
  B[1] = -det*a21;
  B[2] = -det*a12;
  B[3] = det*a11;
  return 1;
}

//=============================================================================
/* Inverse a 3x3 matrix.
   IN: A: 3x3 matrix
   OUT: B: 3x3 matrix
   Retourne 0 si la matrice n'est pas inversible */
//=============================================================================
E_Int K_LINEAR::inv3(E_Float* A, E_Float* B)
{
  E_Float det;
  E_Float a11 = A[0]; E_Float a21 = A[1]; E_Float a31 = A[2]; 
  E_Float a12 = A[3]; E_Float a22 = A[4]; E_Float a32 = A[5];
  E_Float a13 = A[6]; E_Float a23 = A[7]; E_Float a33 = A[8];

  det = a11*(a33*a22-a32*a23)-a21*(a33*a12-a32*a13)+a31*(a23*a12-a22*a13);

  if (K_FUNC::E_abs(det) < 1.e-12) 
  { B[0] = 0.; B[1] = 0.; B[2] = 0.; B[3] = 0.; B[4] = 0.; 
    B[5] = 0.; B[6] = 0.; B[7] = 0.; return 0; }

  det = 1./det;
  B[0] = det*(a33*a22-a32*a23);
  B[1] = -det*(a33*a21-a31*a23);
  B[2] = det*(a32*a21-a31*a22);
  B[3] = -det*(a33*a12-a32*a13);
  B[4] = det*(a33*a11-a31*a13);
  B[5] = -det*(a32*a11-a31*a12);
  B[6] = det*(a23*a12-a22*a13);
  B[7] = -det*(a23*a11-a21*a13);
  B[8] = det*(a22*a11-a21*a12); 
  return 1;
}
