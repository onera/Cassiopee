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

#include "snode.h"

struct Status {
    Snode *root;
    E_Float rx, ry;

    Status();

    Snode *insert(Segment *key, void *inf = NULL);

    Snode *lookup(Segment *seg);

    Snode *locate(Segment *seg);

    Snode *pred(Snode *sit);
    
    Snode *pred(Segment *seg);
    
    Snode *succ(Snode *sit);
    
    Snode *succ(Segment *seg);

    void erase(Segment *seg);

    void print();
    
    Snode *insert_(Snode *& root, Segment *key, void *inf);
    
    Snode *lookup_(Snode *root, Segment *seg);

    Snode *locate_(Snode *root, Segment *seg);

    void pred_(Snode *root, Segment *seg, Snode *&pre);
    
    void succ_(Snode *root, Segment *seg, Snode *&suc);

    Snode *erase_(Snode *root, Segment *seg);
};
