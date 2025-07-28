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
#include "math.h"

const E_Float K_MATH::ONE_THIRD = 0.333333333333333;
const E_Float K_MATH::PI        = 3.141592653589793;
const E_Float K_MATH::SMALL     = 1.0e-15;

void K_MATH::sqrmat_dot_vec(const E_Float *A, const E_Float *x, const E_Int n,
  E_Float *y)
{
  for (E_Int i = 0; i < n; i++)
  {
    y[i] = 0.0;
    const E_Float *ptr = &A[i*n];
    for (E_Int j = 0; j < n; j++) y[i] += ptr[j]*x[j];
  }
}

void K_MATH::sym3mat_dot_vec(const E_Float a[6], const E_Float b[3],
  E_Float c[3])
{
  c[0] = a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
  c[1] = a[1]*b[0] + a[3]*b[1] + a[4]*b[2];
  c[2] = a[2]*b[0] + a[4]*b[1] + a[5]*b[2];
}

E_Float K_MATH::sign(const E_Float a, const E_Float tol)
{
  if (K_MATH::feq(a, 0.0, tol)) return 0.0;
  return a > 0.0 ? 1.0 : -1.0;
}

void K_MATH::sym3mat_dot_sym3mat(const E_Float a[6], const E_Float b[6],
  E_Float c[9])
{
  c[0] = a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
  c[1] = a[0]*b[1] + a[1]*b[3] + a[2]*b[4];
  c[2] = a[0]*b[2] + a[1]*b[4] + a[2]*b[5]; 
  c[3] = a[1]*b[0] + a[3]*b[1] + a[4]*b[2];
  c[4] = a[1]*b[1] + a[3]*b[3] + a[4]*b[4];
  c[5] = a[1]*b[2] + a[3]*b[4] + a[4]*b[5]; 
  c[6] = a[2]*b[0] + a[4]*b[1] + a[5]*b[2];
  c[7] = a[2]*b[1] + a[4]*b[3] + a[5]*b[4];
  c[8] = a[2]*b[2] + a[4]*b[4] + a[5]*b[5];
}

E_Float K_MATH::sym3mat_det(const E_Float A[6])
{
  return A[0]*A[3]*A[5] - (A[0]*A[4]*A[4] + A[3]*A[2]*A[2] + A[5]*A[1]*A[1]) +
    2.0*A[1]*A[2]*A[4];
}

E_Float K_MATH::sym3mat_trace(const E_Float A[6])
{
  return A[0] + A[3] + A[5];
}

E_Float K_MATH::sym3mat_second_invariant(const E_Float A[6])
{
  E_Float AA[9];
  sym3mat_dot_sym3mat(A, A, AA);
  return 0.5*(AA[0] + AA[4] + AA[8]);
}

E_Float K_MATH::sym3mat_third_invariant(const E_Float A[6])
{
  return sym3mat_det(A);
}
