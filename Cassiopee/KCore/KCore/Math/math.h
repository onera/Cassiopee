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
#ifndef _KCORE_MATH_H
#define _KCORE_MATH_H

#include <cmath>
#include "Def/DefTypes.h"

namespace K_MATH
{
  extern const E_Float ONE_THIRD;
  extern const E_Float PI;
  extern const E_Float SMALL;

  // Compute the dot product of two vectors
  inline void dot(
    const E_Float a1, const E_Float a2, const E_Float a3,
    const E_Float b1, const E_Float b2, const E_Float b3,
    E_Float& res
  )
  {
    res = a1 * b1 + a2 * b2 + a3 * b3;
  }

  inline E_Float dot(const E_Float* a, const E_Float* b, const E_Int n)
  {
    E_Float res = 0.0;
    for (E_Int i = 0; i < n; i++) res += a[i] * b[i];
    return res;
  }
  
  // Compute the cross product of two vectors
  inline void cross(
    const E_Float a1, const E_Float a2, const E_Float a3,
    const E_Float b1, const E_Float b2, const E_Float b3,
    E_Float& c1, E_Float& c2, E_Float& c3
  )
  {
    c1 = a2 * b3 - a3 * b2;
    c2 = a3 * b1 - a1 * b3;
    c3 = a1 * b2 - a2 * b1;
  }

  inline void cross(const E_Float a[3], const E_Float b[3], E_Float c[3])
  {
    c[0] = a[1] * b[2] - a[2] * b[1];
    c[1] = a[2] * b[0] - a[0] * b[2];
    c[2] = a[0] * b[1] - a[1] * b[0];
  }

  inline E_Float norm(const E_Float* a, const E_Int n)
  {
    return sqrt(K_MATH::dot(a, a, n));
  }

  inline E_Int feq(const E_Float a, const E_Float b, const E_Float tol = 1e-12)
  {
    return fabs(a - b) < tol;
  }

  void sqrmat_dot_vec(const E_Float*, const E_Float*, const E_Int, E_Float*);
    
  void sym3mat_dot_vec(const E_Float*, const E_Float*, E_Float*);
  
  E_Float sign(const E_Float, const E_Float tol = 1e-15);
  
  void sym3mat_dot_sym3mat(const E_Float*, const E_Float*, E_Float*);
  
  E_Float sym3mat_det(const E_Float [6]);
  
  E_Float sym3mat_trace(const E_Float [6]);
  
  E_Float sym3mat_second_invariant(const E_Float*);
  
  E_Float sym3mat_third_invariant(const E_Float*);
}

#endif
