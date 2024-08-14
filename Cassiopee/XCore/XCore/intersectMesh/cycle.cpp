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
#include <cstddef>
#include <cassert>
#include <cstdio>
#include <set>

#include "cycle.h"
#include "dcel.h"
#include "hedge.h"
#include "primitives.h"

Int Cycle::INNER = 0;
Int Cycle::OUTER = 1;
Int Cycle::DEGEN = 2;

Cycle::Cycle(Hedge *Rep)
: rep(Rep), inout(Dcel::NO_IDEA), left(NULL), prev(NULL), next(NULL)
{}

void Cycle::write_vertices(const char *fname, const std::vector<Cycle *> &C)
{
    FILE *fh = fopen(fname, "w");
    assert(fh);

    std::set<Vertex *> V;
    for (Cycle *c : C) {
        Hedge *h = c->rep;
        Vertex *v = h->orig;
        if (V.find(v) == V.end()) V.insert(v);

        Hedge *w = h->next;
        while (w != h) {
            v = w->orig;
            if (V.find(v) == V.end()) V.insert(v);
            w = w->next;
        }
    }

    fprintf(fh, "POINTS\n");
    fprintf(fh, "%zu\n", V.size());
    for (const auto &v : V) {
        fprintf(fh, "%f %f %f\n", v->x, v->y, v->z);
    }
    

    fclose(fh);
}
