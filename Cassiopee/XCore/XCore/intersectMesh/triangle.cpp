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

E_Int Triangle::is_point_inside(E_Float px, E_Float py, E_Float pz,
    E_Float ax, E_Float ay, E_Float az,
    E_Float bx, E_Float by, E_Float bz,
    E_Float cx, E_Float cy, E_Float cz,
    E_Float &u, E_Float &v, E_Float &w)
{
    // Normal vector to the plane
    E_Float v0[3] = {bx-ax, by-ay, bz-az};
    E_Float v1[3] = {cx-ax, cy-ay, cz-az};
    E_Float v2[3] = {px-ax, py-ay, pz-az};
    
    E_Float N[3];
    K_MATH::cross(v0, v1, N);

    // TODO(Imad): check degenerate triangle with |N| = 0

    E_Float dp = K_MATH::dot(N, v2, 3);
    
    // Is the point on the plane?
    if (dp < -TOL || dp > TOL) return 0;

    E_Float d00 = K_MATH::dot(v0, v0, 3);
    E_Float d01 = K_MATH::dot(v0, v1, 3);
    E_Float d11 = K_MATH::dot(v1, v1, 3);
    E_Float d20 = K_MATH::dot(v2, v0, 3);
    E_Float d21 = K_MATH::dot(v2, v1, 3);

    E_Float denom = d00 * d11 - d01 * d01;

    v = (d11*d20 - d01*d21) / denom;
    if (v < -TOL || v > 1 + TOL) return 0;
    
    w = (d00*d21 - d01*d20) / denom;
    if (w < -TOL || w > 1 + TOL) return 0;

    u = 1 - v - w;
    if (u < -TOL || u > 1 + TOL) return 0;

    return 1;
}

E_Int Triangle::is_point_inside(E_Float px, E_Float py, E_Float pz,
    E_Float ax, E_Float ay, E_Float az,
    E_Float bx, E_Float by, E_Float bz,
    E_Float cx, E_Float cy, E_Float cz)
{
    // Normal vector to the plane
    E_Float v0[3] = {bx-ax, by-ay, bz-az};
    E_Float v1[3] = {cx-ax, cy-ay, cz-az};
    E_Float v2[3] = {px-ax, py-ay, pz-az};
    
    E_Float N[3];
    K_MATH::cross(v0, v1, N);

    // TODO(Imad): check degenerate triangle with |N| = 0

    E_Float dp = K_MATH::dot(N, v2, 3);
    
    // Is the point on the plane?
    if (dp < -TOL || dp > TOL) return 0;

    E_Float d00 = K_MATH::dot(v0, v0, 3);
    E_Float d01 = K_MATH::dot(v0, v1, 3);
    E_Float d11 = K_MATH::dot(v1, v1, 3);
    E_Float d20 = K_MATH::dot(v2, v0, 3);
    E_Float d21 = K_MATH::dot(v2, v1, 3);

    E_Float denom = d00 * d11 - d01 * d01;

    E_Float v = (d11*d20 - d01*d21) / denom;
    if (v < -TOL || v > 1 + TOL) return 0;
    
    E_Float w = (d00*d21 - d01*d20) / denom;
    if (w < -TOL || w > 1 + TOL) return 0;

    E_Float u = 1 - v - w;
    if (u < -TOL || u > 1 + TOL) return 0;

    return 1;
}
