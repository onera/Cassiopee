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
#include <cstdio>
#include <cassert>
#include <cmath>

#include "primitives.h"

E_Int Sign(E_Float x, E_Float tol)
{
    if (x > tol) return 1;
    if (x < -tol) return -1;
    return 0;
}

E_Int cmp_points(E_Float x1, E_Float y1, E_Float z1,
    E_Float x2, E_Float y2, E_Float z2)
{
    if (x1 < x2) return -1;
    if (x1 > x2) return 1;
    if (y1 < y2) return -1;
    if (y1 > y2) return 1;
    if (z1 < z2) return -1;
    if (z1 > z2) return 1;
    return 0;
    /*
    E_Float t = x1 - x2;
    E_Int s = Sign(t);
    if (s) return s;
    
    t = y1 - y2;
    s = Sign(t);
    if (s) return s;

    t = z1 - z2;
    return Sign(t);
    */
}

bool ray_edge_intersect(E_Float ox, E_Float oy, E_Float oz,
    E_Float dx, E_Float dy, E_Float dz,
    E_Float px, E_Float py, E_Float pz,
    E_Float qx, E_Float qy, E_Float qz,
    E_Float &t, E_Float &u,
    bool coplanar)
{
    E_Float v[3]= {px-ox, py-oy, pz-oz};

    E_Float dl[3] = {qx-px, qy-py, qz-pz};

    E_Float dr[3] = {dx, dy, dz};

    if (!coplanar) {
        E_Float w[3];
        K_MATH::cross(v, dl, w);
        E_Float det = K_MATH::dot(w, dr, 3);

        // Ray and edge must be coplanar
        if (Sign(det) != 0) return false;
    }

    // ray and edge must not be parallel
    E_Float n[3];
    K_MATH::cross(dr, dl, n);
    E_Float denom = K_MATH::dot(n, n, 3);
    if (Sign(denom) == 0) return false;

    E_Float tmp[3];
    K_MATH::cross(v, dl, tmp);

    t = K_MATH::dot(tmp, n, 3) / denom;

    if (t < TOL) return false;

    K_MATH::cross(v, dr, tmp);

    u = K_MATH::dot(tmp, n, 3) / denom;

    if (u < -RAY_EDGE_TOL || u > 1 + RAY_EDGE_TOL) return false;

    return true;
}
