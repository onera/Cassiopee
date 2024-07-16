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
#include "common/mem.h"

int Mesh_conformize_cell_face(Mesh *M, Int cid, Int fid, Int fpos, Int nf)
{
    Int *cell = Mesh_get_cell(M, cid);
    Int *crange = Mesh_get_crange(M, cid);

    Int found = 0;

    for (Int i = 0; i < M->cstride[cid]; i++) {
        Int *pf = cell + 4*i;
 
        if (pf[0] == fid) {
            crange[i] = nf;
            for (Int j = 1; j < nf; j++) {
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

    for (Int fid = 0; fid < M->nf; fid++) {
        Int *fpts = Mesh_get_face(M, fid);
        Int *frange = Mesh_get_frange(M, fid);
        Int size = 2*M->fstride[fid];

        for (Int i = 0; i < size; i += 2) {
            Int p = fpts[i];
            Int q = fpts[(i+2)%size];

            UEdge E(p, q);
            
            auto it = M->ecenter.find(E);

            if (it != M->ecenter.end()) {
                Int pid = it->second;
                fpts[i+1] = pid;
                frange[i/2] = 2;
            }
        }
    }

    // Exchange franges at the interface

    for (Int i = 0; i < M->npp; i++) {
        PPatch *P = &M->pps[i];

        P->sbuf_i = (Int *)XRESIZE(P->sbuf_i, 4 * P->nf * sizeof(Int));
        P->rbuf_i = (Int *)XRESIZE(P->rbuf_i, 4 * P->nf * sizeof(Int));

        memset(P->sbuf_i, 0, 4*P->nf * sizeof(Int));
        memset(P->rbuf_i, 0, 4*P->nf * sizeof(Int));

        for (Int j = 0; j < P->nf; j++) {
            Int fid = P->pf[j];
            Int *frange = Mesh_get_frange(M, fid);

            Int *ptr = P->sbuf_i + 4*j;

            for (Int k = 0; k < 4; k++) ptr[k] = frange[k];

            //Int fstate = M->fref[fid];
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

    Int pcount = 0;

    for (Int i = 0; i < M->npp; i++) {
        PPatch *P = &M->pps[i];

        for (Int j = 0; j < P->nf; j++) {
            Int *range = P->rbuf_i + 4*j;

            for (Int k = 0; k < 4; k++) {
                if (range[k] == 2) pcount++;
            }
        }
    }

    Int plimit = M->np + pcount;

    M->X = (Float *)XRESIZE(M->X, plimit * sizeof(Float));
    M->Y = (Float *)XRESIZE(M->Y, plimit * sizeof(Float));
    M->Z = (Float *)XRESIZE(M->Z, plimit * sizeof(Float));

    // Create new edge-center entries if needed

    std::map<UEdge, Int> new_ecenter;

    for (Int i = 0; i < M->npp; i++) {
        PPatch *P = &M->pps[i];

        for (Int j = 0; j < P->nf; j++) {
            Int *range = P->rbuf_i + 4*j;
            Int fid = P->pf[j];
            Int *face = Mesh_get_face(M, fid);
            Int *frange = Mesh_get_frange(M, fid);

            Int size = 2*M->fstride[fid];

            for (Int k = 0; k < size; k += 2) {
                if (range[k/2] == 2 && frange[k/2] == 1) {
                    Int p = face[k];
                    Int q = face[(k+2)%size];

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

    for (Int fid = 0; fid < M->nf; fid++) {
        Int *fpts = Mesh_get_face(M, fid);
        Int *frange = Mesh_get_frange(M, fid);
        Int size = 2*M->fstride[fid];

        for (Int i = 0; i < size; i += 2) {
            if (frange[i/2] == 2) continue;

            Int p = fpts[i];
            Int q = fpts[(i+2)%size];

            UEdge E(p, q);
            
            auto it = new_ecenter.find(E);

            if (it != new_ecenter.end()) {
                Int pid = it->second;
                fpts[i+1] = pid;
                frange[i/2] = 2;
            }
        }
    }
}