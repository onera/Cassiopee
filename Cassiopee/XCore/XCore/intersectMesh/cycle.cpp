#include <cstddef>
#include <cassert>
#include <cstdio>
#include <set>

#include "cycle.h"
#include "dcel.h"
#include "hedge.h"
#include "primitives.h"

E_Int Cycle::INNER = 0;
E_Int Cycle::OUTER = 1;
E_Int Cycle::DEGEN = 2;

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
    fprintf(fh, "%lu\n", V.size());
    for (const auto &v : V) {
        fprintf(fh, "%f %f 0.0\n", v->x, v->y);
    }
    

    fclose(fh);
}
