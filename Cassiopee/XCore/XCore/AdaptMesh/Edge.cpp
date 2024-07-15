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
