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

    E_Int nelem;

    Queue();

    Event *insert(E_Float X, E_Float Y, E_Float Z, E_Int oid, E_Int color);

    Event *insert(E_Float X, E_Float Y, E_Float Z);

    Event *lookup(Vertex *key);
    
    Event *lookup(E_Float x, E_Float y, E_Float z);
    
    inline bool empty() { return root == NULL; }

    Event *min();

    void erase(Event *event);

    void erase(Vertex *p);

    void inorder(std::vector<Vertex *> &V) const;
    
    Event *insert_(Event *&root, E_Float x, E_Float y, E_Float z, E_Int oid, E_Int color);

    Event *insert_(Event *&root, E_Float x, E_Float y, E_Float z);

    Event *lookup_(Event *root, E_Float x, E_Float y, E_Float z);

    Event *erase_(Event *root, Vertex *p);

    void drop();
};
