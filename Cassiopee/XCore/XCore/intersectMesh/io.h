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
#include <unordered_map>

#include "xcore.h"
#include "vertex.h"
#include "triangleIntersection.h"
#include "common/common.h"
#include "point.h"

void point_write(const char *fname, const std::vector<point> &P);

void point_write(const char *fname, const std::vector<Vertex *> &I);

void point_write(const char *fname, Float *Xs, Float *Ys, Float *Zs,
    const std::vector<Int> &proj_points);

void edge_write(const char *fname, Float *X, Float *Y, Float *Z,
    const std::unordered_map<Int, TriangleIntersection> &point_hits);
