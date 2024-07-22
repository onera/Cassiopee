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

#include "xcore.h"
#include "common/common.h"

struct Hedge;
struct Vertex;

struct Cycle {
    Hedge *rep;
    Int inout;
    Vertex *left;
    Cycle *prev;
    Cycle *next;

    static Int INNER;
    static Int OUTER;
    static Int DEGEN;

    Cycle(Hedge *Rep);
    void print() const;

    static void set_inout(std::vector<Cycle *> &C);

    static void write_vertices(const char *fname, const std::vector<Cycle *> &C);
};