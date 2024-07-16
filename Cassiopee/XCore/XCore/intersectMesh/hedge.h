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

struct Vertex;
struct Face;
struct Cycle;

struct Hedge {
    Vertex *orig;
    Hedge *twin;
    Hedge *prev;
    Hedge *next;
    Face *left;
    E_Int color;
    Cycle *cycle;

    Hedge(Vertex *Orig);

    static E_Int cmp_cwise(const Hedge *h, const Hedge *w);
    static void sort_cwise(std::vector<Hedge *> &H, E_Int start, E_Int end);
    static void sort_ccwise(std::vector<Hedge *> &H, E_Int start, E_Int end);
};