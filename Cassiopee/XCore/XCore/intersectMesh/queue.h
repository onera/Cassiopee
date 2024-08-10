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

#include <cstddef>
#include <vector>

#include "xcore.h"
#include "common/common.h"

struct Vertex;
struct Segment;
struct Event;

struct Queue {
    Event *root;

    Int nelem;

    Queue();

    Event *insert(Float X, Float Y, Float Z, Int oid, Int color);

    Event *insert(Float X, Float Y, Float Z);

    Event *lookup(Vertex *key);
    
    Event *lookup(Float x, Float y, Float z);
    
    inline bool empty() { return root == NULL; }

    Event *min();

    void erase(Event *event);

    void erase(Vertex *p);

    void inorder(std::vector<Vertex *> &V) const;
    
    Event *insert_(Event *&root, Float x, Float y, Float z, Int oid, Int color);

    Event *insert_(Event *&root, Float x, Float y, Float z);

    Event *lookup_(Event *root, Float x, Float y, Float z);

    Event *erase_(Event *root, Vertex *p);
};
