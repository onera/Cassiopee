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

//Authors : Sam Landier (sam.landier@onera.fr)

#ifndef __LINEAR_DELAUNAY_MATH_H__
#define __LINEAR_DELAUNAY_MATH_H__

#include "Nuga/include/defs.h"
#include "Nuga/include/DynArray.h"

namespace K_LINEAR
{
  class DelaunayMath
  {
  public:
    static void eigen_vectors 
      (E_Float a00, E_Float a11, E_Float a10, E_Float& lambda0, E_Float& lambda1, E_Float* v0, E_Float* v1);

    static void eigen_values 
      (E_Float a00, E_Float a11, E_Float a10, E_Float& lambda0, E_Float& lambda1);

    static void simultaneous_reduction
      (const K_FLD::FloatArray& M1, const K_FLD::FloatArray& M2,
       E_Float* V1, E_Float* V2);

    static void intersect(const K_FLD::FloatArray& Mi, const K_FLD::FloatArray& M2, K_FLD::FloatArray& I);

    static void resoLin(E_Float a11, E_Float a12, E_Float a21, E_Float a22,
                        E_Float d1, E_Float d2, E_Float&x, E_Float&y);
  };
}

#endif
