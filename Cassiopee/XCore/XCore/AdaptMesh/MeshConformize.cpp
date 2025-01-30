/*    
    Copyright 2013-2025 Onera.

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
#include "common/mem.h"

E_Int Mesh_conformize_cell_face(Mesh *M, E_Int cid, E_Int fid, E_Int fpos, E_Int nf)
{
    E_Int *cell = Mesh_get_cell(M, cid);
    E_Int *crange = Mesh_get_crange(M, cid);

    E_Int found = 0;

    for (E_Int i = 0; i < M->cstride[cid]; i++) {
        E_Int *pf = cell + 4*i;
 
        if (pf[0] == fid) {
            crange[i] = nf;
            for (E_Int j = 1; j < nf; j++) {
                pf[j] = fpos+j-1;
            }
            found = 1;
            break;
        }
    }

    return !found;
}


// TODO(Imad): optimize

void Mesh_conformize_face_edge(Mesh *M)
{
    // TODO(Imad): is it ok to visit only old faces?

    for (E_Int fid = 0; fid < M->nf; fid++) {
        E_Int *fpts = Mesh_get_face(M, fid);
        E_Int *frange = Mesh_get_frange(M, fid);
        E_Int size = 2*M->fstride[fid];

        for (E_Int i = 0; i < size; i += 2) {
            E_Int p = fpts[i];
            E_Int q = fpts[(i+2)%size];

            UEdge E(p, q);
            
            auto it = M->ecenter.find(E);

            if (it != M->ecenter.end()) {
                E_Int pid = it->second;
                fpts[i+1] = pid;
                frange[i/2] = 2;
            }
        }
    }

    if (M->npc == 1) return;

    // Exchange franges at the interface

    for (E_Int i = 0; i < M->npp; i++) {
        PPatch *P = &M->pps[i];

        P->sbuf_i = (E_Int *)XRESIZE(P->sbuf_i, 4 * P->nf * sizeof(E_Int));
        P->rbuf_i = (E_Int *)XRESIZE(P->rbuf_i, 4 * P->nf * sizeof(E_Int));

        memset(P->sbuf_i, 0, 4*P->nf * sizeof(E_Int));
        memset(P->rbuf_i, 0, 4*P->nf * sizeof(E_Int));

        for (E_Int j = 0; j < P->nf; j++) {
            E_Int fid = P->pf[j];
            E_Int *frange = Mesh_get_frange(M, fid);

            E_Int *ptr = P->sbuf_i + 4*j;

            for (E_Int k = 0; k < 4; k++) ptr[k] = frange[k];

            //E_Int fstate = M->fref[fid];
            //assert(fstate != FACE_NEW);

            // Face is oriented in the opposite order
            std::reverse(ptr, ptr+4);

            //k += 4;
        }

        MPI_Isend(P->sbuf_i, 4*P->nf, XMPI_INT, P->nei, M->pid,
            MPI_COMM_WORLD, &M->reqs[M->nrq++]);
        MPI_Irecv(P->rbuf_i, 4*P->nf, XMPI_INT, P->nei, P->nei,
            MPI_COMM_WORLD, &M->reqs[M->nrq++]);
    }

    Mesh_comm_waitall(M);

    // Max estimation of number of new points

    E_Int pcount = 0;

    for (E_Int i = 0; i < M->npp; i++) {
        PPatch *P = &M->pps[i];

        for (E_Int j = 0; j < P->nf; j++) {
            E_Int *range = P->rbuf_i + 4*j;

            for (E_Int k = 0; k < 4; k++) {
                if (range[k] == 2) pcount++;
            }
        }
    }

    E_Int plimit = M->np + pcount;

    M->X = (E_Float *)XRESIZE(M->X, plimit * sizeof(E_Float));
    M->Y = (E_Float *)XRESIZE(M->Y, plimit * sizeof(E_Float));
    M->Z = (E_Float *)XRESIZE(M->Z, plimit * sizeof(E_Float));

    // Create new edge-center entries if needed

    std::map<UEdge, E_Int> new_ecenter;

    for (E_Int i = 0; i < M->npp; i++) {
        PPatch *P = &M->pps[i];

        for (E_Int j = 0; j < P->nf; j++) {
            E_Int *range = P->rbuf_i + 4*j;
            E_Int fid = P->pf[j];
            E_Int *face = Mesh_get_face(M, fid);
            E_Int *frange = Mesh_get_frange(M, fid);

            E_Int size = 2*M->fstride[fid];

            for (E_Int k = 0; k < size; k += 2) {
                if (range[k/2] == 2 && frange[k/2] == 1) {
                    E_Int p = face[k];
                    E_Int q = face[(k+2)%size];

                    UEdge E(p, q);

                    auto it = M->ecenter.find(E);

                    // Edge already conformized, skip
                    if (it != M->ecenter.end()) continue;

                    it = new_ecenter.find(E);
                    
                    if (it == new_ecenter.end()) {
                        assert(M->np < plimit);
                        M->X[M->np] = 0.5 * (M->X[p] + M->X[q]);
                        M->Y[M->np] = 0.5 * (M->Y[p] + M->Y[q]);
                        M->Z[M->np] = 0.5 * (M->Z[p] + M->Z[q]);
                        new_ecenter[E] = M->np++;
                    }
                }
            }
        }
    }

    // Conformize again if needed

    if (new_ecenter.empty()) return;

    for (E_Int fid = 0; fid < M->nf; fid++) {
        E_Int *fpts = Mesh_get_face(M, fid);
        E_Int *frange = Mesh_get_frange(M, fid);
        E_Int size = 2*M->fstride[fid];

        for (E_Int i = 0; i < size; i += 2) {
            if (frange[i/2] == 2) continue;

            E_Int p = fpts[i];
            E_Int q = fpts[(i+2)%size];

            UEdge E(p, q);
            
            auto it = new_ecenter.find(E);

            if (it != new_ecenter.end()) {
                E_Int pid = it->second;
                fpts[i+1] = pid;
                frange[i/2] = 2;
            }
        }
    }
}