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

#include "queue.h"
#include "status.h"
#include "common/common.h"

extern E_Float TOL;

void compute_intersection(Queue &Q, Snode *sit0, Snode *sit1,
    std::vector<Vertex *> &I);

E_Int compare(const Vertex &a, const Vertex &b);

E_Int compare(const Segment &s0, const Segment &s1, E_Float rx, E_Float ry);

E_Int cmp_mySeg(const Segment &s1, const Segment &s2);

E_Int cmp_points(E_Float x1, E_Float y1, E_Float z1, E_Float x2, E_Float y2, E_Float z2);
//E_Int cmp_points(E_Float x1, E_Float y1, E_Float x2, E_Float y2);

E_Float DifferenceOfProducts(E_Float a, E_Float b, E_Float c, E_Float d);

E_Float TwoDiff(E_Float a, E_Float b);

E_Int Sign(E_Float x);

E_Int orient3D(E_Float *A, E_Float *B, E_Float *C, E_Float *D);

E_Float dRand(E_Float dMin, E_Float dMax);

E_Int is_point_on_segment(E_Float px, E_Float py, E_Float pz, E_Float ax, E_Float ay,
    E_Float az, E_Float bx, E_Float by, E_Float bz);

E_Int EdgeEdgeIntersect(E_Float ax, E_Float ay, E_Float az, E_Float bx, E_Float by,
    E_Float bz, E_Float px, E_Float py, E_Float pz, E_Float qx, E_Float qy, E_Float qz,
    E_Float &ix, E_Float &iy, E_Float &iz);

E_Int EdgeEdgeIntersect(E_Float ax, E_Float ay, E_Float az, E_Float bx, E_Float by,
    E_Float bz, E_Float px, E_Float py, E_Float pz, E_Float qx, E_Float qy, E_Float qz);

E_Int EdgeEdgeIntersect(E_Float ax, E_Float ay, E_Float az, E_Float bx, E_Float by,
    E_Float bz, E_Float px, E_Float py, E_Float pz, E_Float qx, E_Float qy, E_Float qz,
    E_Float &t);
