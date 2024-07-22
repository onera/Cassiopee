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

extern Float TOL;

void compute_intersection(Queue &Q, Snode *sit0, Snode *sit1,
    std::vector<Vertex *> &I);

Int compare(const Vertex &a, const Vertex &b);

Int compare(const Segment &s0, const Segment &s1, Float rx, Float ry);

Int cmp_mySeg(const Segment &s1, const Segment &s2);

Int cmp_points(Float x1, Float y1, Float x2, Float y2);

Float DifferenceOfProducts(Float a, Float b, Float c, Float d);

Float TwoDiff(Float a, Float b);

Int Sign(Float x);

Int orient3D(Float *A, Float *B, Float *C, Float *D);

Float dRand(Float dMin, Float dMax);
