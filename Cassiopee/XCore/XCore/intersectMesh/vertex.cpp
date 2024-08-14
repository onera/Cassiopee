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
#include "vertex.h"
#include "primitives.h"

Vertex::Vertex(Float X, Float Y, Float Z, Int Oid, Int color)
: x(X), y(Y), z(Z), rep(NULL), id(-1), left(NULL)
{
    oid[color] = Oid;
    oid[(color+1)%2] = -1;
}

Vertex::Vertex(Float X, Float Y, Float Z)
: x(X), y(Y), z(Z), rep(NULL), id(-1), left(NULL)
{
    oid[0] = -1;
    oid[1] = -1;
}

Int cmp_vtx(Vertex *a, Vertex *b)
{
    return cmp_points(a->x, a->y, a->z, b->x, b->y, b->z);
}
