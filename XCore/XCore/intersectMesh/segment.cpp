#include "proto.h"
#include <algorithm>
#include <cassert>

segment::segment(hedge *h, int i)
: p(h->orig), q(h->twin->orig), id(i), rep(h)
{
    int cmp = vertex_cmp(p, q);

    if (cmp > 0) {
        vertex *tmp = p;
        p = q;
        q = tmp;
        rep = h->twin;
    }

    assert(sign(q->x - p->x) >= 0);
}

segment::segment(edge *e)
: p(e->p), q(e->q), id(-1), rep(NULL)
{
    int cmp = vertex_cmp(p, q);

    if (cmp > 0) {
        vertex *tmp = p;
        p = q;
        q = tmp;
    }
}

segment::segment(vertex *P, vertex *Q, int Id)
: p(P), q(Q), id(Id), rep(NULL)
{
    int cmp = vertex_cmp(p, q);

    if (cmp > 0) {
        vertex *tmp = p;
        p = q;
        q = tmp;
    }
}

int segment_cmp_lexico(segment *s0, segment *s1)
{
    int cmp = vertex_cmp(s0->p, s1->p);
    if (cmp)
        return cmp;
    cmp = sign(s0->color() - s1->color());
    if (cmp)
        return cmp;
    return sign(s0->id - s1->id);
}

int _partition(std::vector<segment *> &S, int (*cmp)(segment *, segment *),
    int low, int high)
{
    segment *pivot = S[high];
    
    int i = low-1;

    for (int j = low; j < high; j++) {
        if (cmp(S[j], pivot) <= 0) {
            i++;
            std::swap(S[i], S[j]);
        }
    }

    i++;
    std::swap(S[i], S[high]);
    return i;
}
        

void segment_sort(std::vector<segment *> &S, int (*cmp)(segment *, segment *),
    int low, int high)
{
    if (low >= high)
        return;

    int part = _partition(S, cmp, low, high);

    segment_sort(S, cmp, low, part-1);
    segment_sort(S, cmp, part+1, high);
}

int segments_are_colli(segment *s0, segment *s1)
{
    return vertex_orient(s0->p, s0->q, s1->p) == 0 &&
           vertex_orient(s0->p, s0->q, s1->q) == 0;
}


