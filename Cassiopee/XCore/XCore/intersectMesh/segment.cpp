#include <algorithm>
#include <cassert>

#include "segment.h"
#include "primitives.h"
#include "hedge.h"
#include "dcel.h"

Segment::Segment(Vertex *P, Vertex *Q, E_Int Id)
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

Segment::Segment(Hedge *h, E_Int Id)
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
    E_Float sdx = s.dx;
    E_Float sdy = s.dy;
    E_Float sqx = s.q->x;
    E_Float sqy = s.q->y;
    E_Float px = p->x;
    E_Float py = p->y;

    E_Float T1 = dy * sdx - sdy * dx;
    E_Int sign1 = Sign(T1);

    if (sign1 == 0) {
        E_Float mdx = sqx - px;
        E_Float mdy = sqy - py;

        E_Int sign2 = Sign(dy * mdx - mdy * dx);

        if (sign2 == 0) {
            E_Int sign3 = Sign(sdy * mdx - mdy * sdx);

            assert(sign3 == 0);

            if (sign3 == 0) {
                return true;
            }
        }
    }

    return false;
}

static
E_Int _partition(std::vector<Segment *> &S, E_Int low, E_Int high,
    E_Int (*cmp)(const Segment &, const Segment &))
{
    Segment *pivot = S[high];
    E_Int i = low-1;

    for (E_Int j = low; j < high; j++) {
        if (cmp(*S[j], *pivot) <= 0) {
            i++;
            std::swap(S[i], S[j]);
        }
    }

    i++;
    std::swap(S[i], S[high]);
    return i;
}

void Segment::sort(std::vector<Segment *> &S, E_Int start, E_Int end,
    E_Int (*cmp)(const Segment &, const Segment &))
{
    if (start >= end) return;

    E_Int p = _partition(S, start, end, cmp);

    sort(S, start, p - 1, cmp);
    sort(S, p + 1, end, cmp);
}
