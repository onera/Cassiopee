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
#pragma once

#include "triangle.h"
#include "point.h"

struct AABB;

struct Vec3 {
    E_Float x, y, z;
};

struct Ray {
    Point org;
    Vec3 dir;
    E_Int kx, ky, kz;
    E_Float Sx, Sy, Sz;

    Ray(Point O, Vec3 D);

    E_Int intersect_triangle(const Point &a, const Point &b, const Point &c,
        TriangleIntersection &TI);
};

E_Int ray_point_orient(const E_Float o[3], const E_Float d[3],
    const E_Float fN[3], E_Float px, E_Float py, E_Float pz);

bool ray_AABB_intersect(E_Float ox, E_Float oy, E_Float oz,
    E_Float dx, E_Float dy, E_Float dz,
    const AABB &box);

E_Int MollerTrumbore(E_Float px, E_Float py, E_Float pz, E_Float dx, E_Float dy,
    E_Float dz, E_Float ax, E_Float ay, E_Float az, E_Float bx, E_Float by, E_Float bz,
    E_Float cx, E_Float cy, E_Float cz, TriangleIntersection &TI);

bool MollerTrumboreAnyDir(
    E_Float px, E_Float py, E_Float pz,
    E_Float dx, E_Float dy, E_Float dz,
    E_Float ax, E_Float ay, E_Float az,
    E_Float bx, E_Float by, E_Float bz,
    E_Float cx, E_Float cy, E_Float cz,
    E_Float &u, E_Float &v, E_Float &w, E_Float &t,
    E_Float &x, E_Float &y, E_Float &z);
