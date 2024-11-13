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

#include "common/common.h"

constexpr E_Float TOL = 1e-11;
constexpr E_Float RAY_EDGE_TOL = 1e-5;

E_Int Sign(E_Float x, E_Float tol=TOL);

E_Int cmp_points(E_Float x1, E_Float y1, E_Float z1,
    E_Float x2, E_Float y2, E_Float z2);

bool ray_edge_intersect(E_Float ox, E_Float oy, E_Float oz,
    E_Float dx, E_Float dy, E_Float dz,
    E_Float px, E_Float py, E_Float pz,
    E_Float qx, E_Float qy, E_Float qz,
    E_Float &t, E_Float &u,
    bool coplanar = true);

E_Int ray_point_orient(const E_Float o[3], const E_Float d[3],
    const E_Float fN[3], E_Float px, E_Float py, E_Float pz);
