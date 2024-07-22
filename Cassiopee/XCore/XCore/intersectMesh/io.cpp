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
#include <cstdio>
#include <cassert>

#include "io.h"

void point_write(const char *fname, const std::vector<Vertex *> &I)
{
    FILE *fh = fopen(fname, "w");
    assert(fh);
    fprintf(fh, "POINTS\n");
    fprintf(fh, "%zu\n", I.size());
    for (auto &v : I) fprintf(fh, "%f %f 0\n", v->x, v->y);
    fclose(fh);
}

void point_write(const char *fname, Float *Xs, Float *Ys, Float *Zs,
    const std::vector<Int> &P)
{
    FILE *fh = fopen(fname, "w");
    assert(fh);
    fprintf(fh, "POINTS\n");
    fprintf(fh, "%zu\n", P.size());
    for (Int p : P) fprintf(fh, "%f %f %f\n", Xs[p], Ys[p], Zs[p]);
    fclose(fh);
}

void edge_write(const char *fname, Float *X, Float *Y, Float *Z,
    const std::unordered_map<Int, TriangleIntersection> &point_hits)
{
    FILE *fh = fopen(fname, "w");
    assert(fh);
    fprintf(fh, "POINTS\n");
    fprintf(fh, "%zu\n", 2 * point_hits.size());
    for (const auto &pdata : point_hits) {
        Int pt = pdata.first;
        const auto &TI = pdata.second;
        fprintf(fh, "%f %f %f\n", X[pt], Y[pt], Z[pt]);
        fprintf(fh, "%f %f %f\n", TI.x, TI.y, TI.z);
    }
    fprintf(fh, "EDGES\n");
    fprintf(fh, "%zu\n", point_hits.size());
    for (size_t i = 0; i < point_hits.size(); i += 2) {
        fprintf(fh, "%zu %zu ", i, i+1);
    }
    fprintf(fh, "\n");
    fclose(fh);
}