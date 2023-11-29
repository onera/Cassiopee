/*    
    Copyright 2013-2024 Onera.

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

#include "Def/DefTypes.h"

namespace K_MATH
{
  extern const E_Float ONE_THIRD;

  extern const E_Float PI;

  extern const E_Float SMALL;

  E_Float dot(const E_Float *, const E_Float *, const E_Int);
  
  void cross(const E_Float [3], const E_Float [3], E_Float [3]);
  
  E_Float norm(const E_Float *, const E_Int);
  
  E_Int feq(const E_Float, const E_Float, const E_Float tol = SMALL);
  
  void sqrmat_dot_vec(const E_Float *, const E_Float *, const E_Int,
    E_Float *);
    
  void sym3mat_dot_vec(const E_Float *, const E_Float *, E_Float *);
  
  E_Float sign(const E_Float, const E_Float tol = 1e-15);
  
  void sym3mat_dot_sym3mat(const E_Float *, const E_Float *, E_Float *);
  
  E_Float sym3mat_det(const E_Float [6]);
  
  E_Float sym3mat_trace(const E_Float [6]);
  
  E_Float sym3mat_second_invariant(const E_Float *);
  
  E_Float sym3mat_third_invariant(const E_Float *);
}

#endif
