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
Int orient(Float ax, Float ay, Float bx, Float by, Float cx, Float cy)
{
    Float det = (bx - ax)*(cy - ay) - (by - ay)*(cx - ax);
    return Sign(det);
}

Int Triangle::is_point_inside(Float px, Float py, Float pz,
    Float ax, Float ay, Float az,
    Float bx, Float by, Float bz,
    Float cx, Float cy, Float cz,
    Float &u, Float &v, Float &w)
{
    // Normal vector to the plane
    Float Y[3] = {bx-ax, by-ay, bz-az};
    Float Z[3] = {cx-ax, cy-ay, cz-az};
    Float N[3];
    K_MATH::cross(Y, Z, N);

    Float X[3] = {px-ax, py-ay, pz-az};

    Float dp = K_MATH::dot(N, X, 3);
    
    // Is the point on the plane?
    if (dp < -TOL || dp > TOL) return 0;

    Float x1 = K_MATH::dot(X, Y, 3);
    Float y1 = K_MATH::dot(Y, Y, 3);
    Float z1 = K_MATH::dot(Z, Y, 3);
    Float x2 = K_MATH::dot(X, Z, 3);
    Float y2 = K_MATH::dot(Y, Z, 3);
    Float z2 = K_MATH::dot(Z, Z, 3);

    u = (x1*z2 - x2*z1) / (y1*z2 - y2*z1);
    if (u < -TOL || u > 1 + TOL) return 0;

    v = (-x1*y2 + x2*y1) / (y1*z2 - y2*z1);
    if (v < -TOL || v > 1 + TOL) return 0;

    w = 1 - u - v;
    if (w < -TOL || w > 1 + TOL) return 0;

    return 1;
}

Int Triangle::is_point_inside(Float px, Float py, Float pz,
    Float ax, Float ay, Float az,
    Float bx, Float by, Float bz,
    Float cx, Float cy, Float cz)
{
    // Normal vector to the plane
    Float Y[3] = {bx-ax, by-ay, bz-az};
    Float Z[3] = {cx-ax, cy-ay, cz-az};
    Float N[3];
    K_MATH::cross(Y, Z, N);

    Float X[3] = {px-ax, py-ay, pz-az};

    Float dp = K_MATH::dot(N, X, 3);
    
    // Is the point on the plane?
    if (dp < -TOL || dp > TOL) return 0;

    Float x1 = K_MATH::dot(X, Y, 3);
    Float y1 = K_MATH::dot(Y, Y, 3);
    Float z1 = K_MATH::dot(Z, Y, 3);
    Float x2 = K_MATH::dot(X, Z, 3);
    Float y2 = K_MATH::dot(Y, Z, 3);
    Float z2 = K_MATH::dot(Z, Z, 3);

    Float u = (x1*z2 - x2*z1) / (y1*z2 - y2*z1);
    if (u < -TOL || u > 1 + TOL) return 0;

    Float v = (-x1*y2 + x2*y1) / (y1*z2 - y2*z1);
    if (v < -TOL || v > 1 + TOL) return 0;

    Float w = 1 - u - v;
    if (w < -TOL || w > 1 + TOL) return 0;

    return 1;
}
