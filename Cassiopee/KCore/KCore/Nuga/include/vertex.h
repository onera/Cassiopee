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

#ifndef NUGA_VERTEX_H
#define NUGA_VERTEX_H

#include "Nuga/include/defs.h"

namespace NUGA
{
  struct vecval
  {
    double  vec[3];
    double  val2;
    int     flag; // e.g. used to store node id
    int     flag2;

    vecval()
    {
      vec[0] = vec[1] = vec[2] = val2 = FLOAT_MAX;
      flag = flag2 = -1;
    }

    vecval(const double*p, double v2):val2(v2)
    {
      vec[0] = p[0];
      vec[1] = p[1];
      vec[2] = p[2];
    }
  };

  using vertex = vecval;
  using direction = vecval;

  
}

#endif
