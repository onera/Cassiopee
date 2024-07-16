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
#include "Mesh.h"

void refine_edge(Int eid, Mesh *M)
{
    Int p = M->edges[eid].p;
    Int q = M->edges[eid].q;
    M->X[M->np] = 0.5 * (M->X[p] + M->X[q]);
    M->Y[M->np] = 0.5 * (M->Y[p] + M->Y[q]);
    M->Z[M->np] = 0.5 * (M->Z[p] + M->Z[q]);

    // Conformize parent faces

    UEdge e(p, q);

    for (Int i = M->xedf[eid]; i < M->xedf[eid+1]; i++) {
        Int fid = M->E2F[i];

        Int *fpts = Mesh_get_face(M, fid);
        Int *frange = Mesh_get_frange(M, fid);

        Int size = 2 * M->fstride[fid];

        Int found = 0;

        for (Int j = 0; j < size; j += 2) {
            Int P = fpts[j];
            Int Q = fpts[(j+2)%size];
            UEdge E(P, Q);
            if (E == e) {
                //assert(fpts[j+1] == -1);
                fpts[j+1] = M->np;
                found = 1;

                frange[j/2] = 2;
                break;
            }
        }

        assert(found);
    }

    M->edges[eid].q = M->np;

    M->np++;
}
