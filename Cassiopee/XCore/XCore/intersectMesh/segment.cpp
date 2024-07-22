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
#include <algorithm>
#include <cassert>

#include "segment.h"
#include "primitives.h"
#include "hedge.h"
#include "dcel.h"

Segment::Segment(Vertex *P, Vertex *Q, Int Id)
{
    p = P;
    q = Q;
    dx = q->x - p->x;
    dy = q->y - p->y;
    rep = NULL;
    id = Id;
    color = Dcel::NO_IDEA;
    assert(Sign(dx) >= 0);
}

Segment::Segment(Vertex *P)
: p(P), q(P), dx(0.0), dy(0.0), rep(NULL), id(-1), color(Dcel::NO_IDEA)
{}

Segment::Segment(Hedge *h, Int Id)
{
    p = h->orig;
    Hedge *t = h->twin;
    q = t->orig;
    rep = h;
    if (compare(*p, *q) > 0) {
        std::swap(p, q);
        rep = t;
    }
    dx = q->x - p->x;
    dy = q->y - p->y;
    id = Id;
    color = rep->color;
    assert(Sign(dx) >= 0);
}

bool Segment::overlaps(const Segment &s)
{
    Float sdx = s.dx;
    Float sdy = s.dy;
    Float sqx = s.q->x;
    Float sqy = s.q->y;
    Float px = p->x;
    Float py = p->y;

    Float T1 = dy * sdx - sdy * dx;
    Int sign1 = Sign(T1);

    if (sign1 == 0) {
        Float mdx = sqx - px;
        Float mdy = sqy - py;

        Int sign2 = Sign(dy * mdx - mdy * dx);

        if (sign2 == 0) {
            Int sign3 = Sign(sdy * mdx - mdy * sdx);

            assert(sign3 == 0);

            if (sign3 == 0) {
                return true;
            }
        }
    }

    return false;
}

static
Int _partition(std::vector<Segment *> &S, Int low, Int high,
    Int (*cmp)(const Segment &, const Segment &))
{
    Segment *pivot = S[high];
    Int i = low-1;

    for (Int j = low; j < high; j++) {
        if (cmp(*S[j], *pivot) <= 0) {
            i++;
            std::swap(S[i], S[j]);
        }
    }

    i++;
    std::swap(S[i], S[high]);
    return i;
}

void Segment::sort(std::vector<Segment *> &S, Int start, Int end,
    Int (*cmp)(const Segment &, const Segment &))
{
    if (start >= end) return;

    Int p = _partition(S, start, end, cmp);

    sort(S, start, p - 1, cmp);
    sort(S, p + 1, end, cmp);
}
