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
#include "Mesh.h"

E_Int refine_edge(E_Int eid, Mesh *M)
{
    E_Int p = M->edges[eid].p;
    E_Int q = M->edges[eid].q;
    M->X[M->np] = 0.5 * (M->X[p] + M->X[q]);
    M->Y[M->np] = 0.5 * (M->Y[p] + M->Y[q]);
    M->Z[M->np] = 0.5 * (M->Z[p] + M->Z[q]);

    // Conformize parent faces

    UEdge e(p, q);

    for (E_Int i = M->xedf[eid]; i < M->xedf[eid+1]; i++) {
        E_Int fid = M->E2F[i];

        E_Int *fpts = Mesh_get_face(M, fid);
        E_Int *frange = Mesh_get_frange(M, fid);

        E_Int size = 2 * M->fstride[fid];

        E_Int found = 0;

        for (E_Int j = 0; j < size; j += 2) {
            E_Int P = fpts[j];
            E_Int Q = fpts[(j+2)%size];
            UEdge E(P, Q);
            if (E == e) {
                //assert(fpts[j+1] == -1);
                fpts[j+1] = M->np;
                found = 1;

                frange[j/2] = 2;
                break;
            }
        }

        if (!found) return 1;
    }

    M->edges[eid].q = M->np;

    M->np++;

    return 0;
}
