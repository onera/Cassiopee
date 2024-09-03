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
#include "hedge.h"

struct IO_Edge {
    E_Float px, py, pz;
    E_Float qx, qy, qz;

    IO_Edge(E_Float PX, E_Float PY, E_Float PZ, E_Float QX, E_Float QY, E_Float QZ)
    : px(PX), py(PY), pz(PZ), qx(QX), qy(QY), qz(QZ)
    {}
};

void face_write(const char *fname, Face *face);

void point_write(const char *fname, E_Float x, E_Float y, E_Float z);

void point_write(const char *fname, Vertex *v);

void hedge_write(const char *fname, const Hedge *h);

void point_write(const char *fname, const std::vector<Point> &P);

void point_write(const char *fname, const std::vector<Vertex *> &I);

void point_write(const char *fname, E_Float *Xs, E_Float *Ys, E_Float *Zs,
    const std::vector<E_Int> &proj_points);

void edge_write(const char *fname, E_Float *X, E_Float *Y, E_Float *Z,
    const std::unordered_map<E_Int, TriangleIntersection> &point_hits);

void edge_write(const char *fname, E_Float px, E_Float py, E_Float pz,
    E_Float qx, E_Float qy, E_Float qz);

void edges_write(const char *fname, const std::vector<IO_Edge> &edges);
