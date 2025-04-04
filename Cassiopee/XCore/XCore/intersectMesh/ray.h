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
#pragma once

#include "triangle.h"
#include "point.h"

struct AABB;

struct Vec3 {
    E_Float x, y, z;
};

struct Ray {
    enum Policy {
        FORWARD = 0, BOTH,
    };

    E_Float o[3];
    E_Float d[3];
    E_Int kx, ky, kz;
    E_Float Sx, Sy, Sz;
    E_Float inv_dx, inv_dy, inv_dz;
    Policy policy;

    Ray(E_Float px, E_Float py, E_Float pz, E_Float dx, E_Float dy, E_Float dz,
        Policy policy);

    bool intersect_AABB(const AABB &box) const;

    bool intersect_triangle(const Point &a, const Point &b, const Point &c,
        E_Float &u, E_Float &v, E_Float &w, E_Float &t, E_Float &x,
        E_Float &y, E_Float &z) const;
};

E_Int ray_point_orient(const E_Float o[3], const E_Float d[3],
    const E_Float fN[3], E_Float px, E_Float py, E_Float pz);

