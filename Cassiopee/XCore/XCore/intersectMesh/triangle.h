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

#include <vector>

#include "point.h"
#include "smesh.h"

struct Triangle {
    Int a, b, c;

    static Int ispointInside(Float x, Float y, Float ax, Float ay,
        Float bx, Float by, Float cx, Float cy);

    static Int ray_intersect(Float px, Float py, Float pz, Float dx,
        Float dy, Float dz, Float ax, Float ay, Float az, Float bx,
        Float by, Float bz, Float cx, Float cy, Float cz, Float &u,
        Float &v, Float &w, Float &t, Float &x, Float &y, Float &z);
};

struct TriIMesh {
    std::vector<point> P;
    std::vector<Triangle> T;

    TriIMesh(const Smesh &M);

    pointFace locate(point p);
};
