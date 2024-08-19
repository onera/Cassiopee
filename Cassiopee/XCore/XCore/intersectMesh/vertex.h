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

#include <cstdio>

#include "xcore.h"
#include "common/common.h"

struct PointLoc {
    E_Int fid = -1;
    E_Int v_idx = -1;
    E_Int e_idx = -1;
};

struct Hedge;

struct Vertex {
    E_Float x, y, z;
    Hedge *rep;
    E_Int id;
    Hedge *left;
    E_Int oid[2];

    PointLoc loc;

    Hedge *xhedge = NULL;

    Vertex(E_Float X, E_Float Y, E_Float Z, E_Int Oid, E_Int color);

    Vertex(E_Float X, E_Float Y, E_Float Z);

    inline void print() { printf("P" SF_D_ ": %f %f %f\n", id, x, y, z); }
};

E_Int cmp_vtx(Vertex *a, Vertex *b);
