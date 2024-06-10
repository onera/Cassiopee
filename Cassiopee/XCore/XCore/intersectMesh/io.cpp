#include <cstdio>
#include <cassert>

#include "io.h"

void point_write(const char *fname, const std::vector<Vertex *> &I)
{
    FILE *fh = fopen(fname, "w");
    assert(fh);
    fprintf(fh, "POINTS\n");
    fprintf(fh, "%lu\n", I.size());
    for (auto &v : I) fprintf(fh, "%f %f 0\n", v->x, v->y);
    fclose(fh);
}

void point_write(const char *fname, E_Float *Xs, E_Float *Ys, E_Float *Zs,
    const std::vector<E_Int> &P)
{
    FILE *fh = fopen(fname, "w");
    assert(fh);
    fprintf(fh, "POINTS\n");
    fprintf(fh, "%lu\n", P.size());
    for (E_Int p : P) fprintf(fh, "%f %f %f\n", Xs[p], Ys[p], Zs[p]);
    fclose(fh);
}