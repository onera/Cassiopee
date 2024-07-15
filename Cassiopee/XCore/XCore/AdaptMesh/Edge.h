#pragma once

#include "../common/common.h" 

struct Mesh;

struct DEdge {
    Int p, q;

    DEdge(Int P, Int Q)
    : p(P), q(Q)
    {}

    bool operator<(const DEdge &E) const
    {
        return (p < E.p) || (p == E.p && q < E.q);
    }
};

void refine_edge(Int eid, Mesh *M);
