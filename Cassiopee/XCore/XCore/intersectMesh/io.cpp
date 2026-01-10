/*    
    Copyright 2013-2026 ONERA.

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

void point_write(const char *fname, E_Float x, E_Float y, E_Float z)
{
    FILE *fh = fopen(fname, "w");
    assert(fh);
    fprintf(fh, "POINTS\n");
    fprintf(fh, "1\n");
    fprintf(fh, "%f %f %f\n", x, y, z);
    fclose(fh);
}

void point_write(const char *fname, E_Float *Xs, E_Float *Ys, E_Float *Zs,
    const std::vector<E_Int> &P)
{
    FILE *fh = fopen(fname, "w");
    assert(fh);
    fprintf(fh, "POINTS\n");
    fprintf(fh, "%zu\n", P.size());
    for (E_Int p : P) fprintf(fh, "%f %f %f\n", Xs[p], Ys[p], Zs[p]);
    fclose(fh);
}

void point_write(const char *fname, const std::vector<Point> &P)
{
    FILE *fh = fopen(fname, "w");
    assert(fh);
    fprintf(fh, "POINTS\n");
    fprintf(fh, "%zu\n", P.size());
    for (auto p : P) fprintf(fh, "%f %f %f\n", p.x, p.y, p.z);
    fclose(fh);
}

void edge_write(const char *fname, E_Float *X, E_Float *Y, E_Float *Z,
    const std::unordered_map<E_Int, TriangleIntersection> &point_hits)
{
    FILE *fh = fopen(fname, "w");
    assert(fh);
    fprintf(fh, "POINTS\n");
    fprintf(fh, "%zu\n", 2 * point_hits.size());
    for (const auto &pdata : point_hits) {
        E_Int pt = pdata.first;
        const auto &TI = pdata.second;
        fprintf(fh, "%f %f %f\n", X[pt], Y[pt], Z[pt]);
        fprintf(fh, "%f %f %f\n", TI.x, TI.y, TI.z);
    }
    fprintf(fh, "EDGES\n");
    fprintf(fh, "%zu\n", point_hits.size());
    for (size_t i = 0; i < 2*point_hits.size(); i++) {
        fprintf(fh, "%zu ", i);
    }
    fprintf(fh, "\n");
    fclose(fh);
}

void edge_write(const char *fname, E_Float px, E_Float py, E_Float pz,
    E_Float qx, E_Float qy, E_Float qz)
{
    FILE *fh = fopen(fname, "w");
    assert(fh);
    fprintf(fh, "POINTS\n");
    fprintf(fh, "2\n");
    fprintf(fh, "%f %f %f\n", px, py, pz);
    fprintf(fh, "%f %f %f\n", qx, qy, qz);
    fprintf(fh, "EDGES\n");
    fprintf(fh, "1\n");
    fprintf(fh, "0 1\n");
    fclose(fh);
}
