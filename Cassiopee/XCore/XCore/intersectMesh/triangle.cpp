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
#include <limits>
#include <cmath>

#include "triangle.h"
#include "primitives.h"

static
E_Int orient(E_Float ax, E_Float ay, E_Float bx, E_Float by, E_Float cx, E_Float cy)
{
    E_Float det = (bx - ax)*(cy - ay) - (by - ay)*(cx - ax);
    return Sign(det);
}

E_Int Triangle::ispointInside(E_Float px, E_Float py, E_Float ax, E_Float ay,
    E_Float bx, E_Float by, E_Float cx, E_Float cy)
{
    if (orient(ax, ay, bx, by, px, py) < 0) return 0;
    if (orient(bx, by, cx, cy, px, py) < 0) return 0;
    if (orient(cx, cy, ax, ay, px, py) < 0) return 0;
    
    return 1;
}