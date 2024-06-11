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

void edge_write(const char *fname, E_Float *X, E_Float *Y, E_Float *Z,
    const std::unordered_map<E_Int, TriangleIntersection> &point_hits)
{
    FILE *fh = fopen(fname, "w");
    assert(fh);
    fprintf(fh, "POINTS\n");
    fprintf(fh, "%lu\n", 2 * point_hits.size());
    for (const auto &pdata : point_hits) {
        E_Int pt = pdata.first;
        const auto &TI = pdata.second;
        fprintf(fh, "%f %f %f\n", X[pt], Y[pt], Z[pt]);
        fprintf(fh, "%f %f %f\n", TI.x, TI.y, TI.z);
    }
    fprintf(fh, "EDGES\n");
    fprintf(fh, "%lu\n", point_hits.size());
    for (size_t i = 0; i < point_hits.size(); i += 2) {
        fprintf(fh, "%lu %lu ", i, i+1);
    }
    fprintf(fh, "\n");
    fclose(fh);
}