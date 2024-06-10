#pragma once

#include <vector>

#include "vertex.h"

struct Hedge;

struct Segment {
    Vertex *p;
    Vertex *q;
    E_Float dx;
    E_Float dy;
    Hedge *rep;
    E_Int id;
    E_Int color;

    inline E_Float X1() const { return p->x; }
    inline E_Float Y1() const { return p->y; }
    inline E_Float X2() const { return q->x; }
    inline E_Float Y2() const { return q->y; }
    
    Segment(Vertex *p, Vertex *q, E_Int Id);
    Segment(Hedge *h, E_Int Id);
    Segment(Vertex *P);

    bool overlaps(const Segment &seg);

    static void sort(std::vector<Segment *> &S, E_Int start, E_Int end,
        E_Int (*cmp)(const Segment &, const Segment &));
};
