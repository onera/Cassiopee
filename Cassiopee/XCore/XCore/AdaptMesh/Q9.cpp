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
#include "Quad.h"
#include "Mesh.h"

E_Int Q9_refine(E_Int quad, Mesh *M)
{
    E_Int NODES[8];

    E_Int *fpts = Mesh_get_face(M, quad);
    E_Int *frange = Mesh_get_frange(M, quad);

    // BOT
    NODES[0] = fpts[0];
    if (frange[0] == 2) {
        NODES[4] = fpts[1];
    } else {
        E_Int p = fpts[0];
        E_Int q = fpts[2];
        Mesh_refine_or_get_edge_center(M, p, q, NODES[4]);
    }

    // RGT
    NODES[1] = fpts[2];
    if (frange[1] == 2) {
        NODES[5] = fpts[3];
    } else {
        E_Int p = fpts[2];
        E_Int q = fpts[4];
        Mesh_refine_or_get_edge_center(M, p, q, NODES[5]);
    }

    // TOP
    NODES[2] = fpts[4];
    if (frange[2] == 2) {
        NODES[6] = fpts[5];
    } else {
        E_Int p = fpts[4];
        E_Int q = fpts[6];
        Mesh_refine_or_get_edge_center(M, p, q, NODES[6]);
    }

    // LFT
    NODES[3] = fpts[6];
    if (frange[3] == 2) {
        NODES[7] = fpts[7];
    } else {
        E_Int p = fpts[6];
        E_Int q = fpts[0];
        Mesh_refine_or_get_edge_center(M, p, q, NODES[7]);
    }

    // Make centroid
    E_Float fc[3] = {0, 0, 0};

    /*
    // Flip to make computation consistent (without FP errors)
    {
        auto it = M->face_to_ppatch.find(quad);
        if (it != M->face_to_bpatch.end()) {
            PPatch *P = &M->pps[it->second];
            if (P->nei < M->pid) {
                E_Int nodes[4];
                for (int i = 0; i < 4; i++) nodes[i] = NODES[i];
                std::reverse(nodes+1, nodes+4);
                for (E_Int i = 0; i < 4; i++) {
                    fc[0] += M->X[nodes[i]];
                    fc[1] += M->Y[nodes[i]];
                    fc[2] += M->Z[nodes[i]];
                }
            }
        }
        else {
            for (E_Int i = 0; i < 4; i++) {
                fc[0] += M->X[NODES[i]];
                fc[1] += M->Y[NODES[i]];
                fc[2] += M->Z[NODES[i]];
            }
        }
    }
    */

    for (E_Int i = 0; i < 4; i++) {
        fc[0] += M->X[NODES[i]];
        fc[1] += M->Y[NODES[i]];
        fc[2] += M->Z[NODES[i]];
    }

    for (E_Int i = 0; i < 3; i++) fc[i] *= 0.25;

    M->X[M->np] = fc[0];
    M->Y[M->np] = fc[1];
    M->Z[M->np] = fc[2];

    E_Int nc = M->np;
    
    // Setup face points

    // Second child
    fpts = Mesh_get_face(M, M->nf);
    fpts[0] = NODES[4]; fpts[2] = NODES[1];
    fpts[4] = NODES[5]; fpts[6] = nc;

    // Third child
    fpts = Mesh_get_face(M, M->nf+1);
    fpts[0] = nc;       fpts[2] = NODES[5];
    fpts[4] = NODES[2]; fpts[6] = NODES[6];

    // Fourth child
    fpts = Mesh_get_face(M, M->nf+2);
    fpts[0] = NODES[7]; fpts[2] = nc;
    fpts[4] = NODES[6]; fpts[6] = NODES[3];

    // First child replaces quad
    fpts = Mesh_get_face(M, quad);
    memset(fpts, -1, 8*sizeof(E_Int));
    fpts[0] = NODES[0]; fpts[2] = NODES[4];
    fpts[4] = nc;       fpts[6] = NODES[7];

    // Update ranges and strides
    Mesh_update_face_range_and_stride(M, quad, M->nf, 3);

    // Conformize parent cells

    E_Int own = M->owner[quad];

    assert(M->clevel[own] == M->flevel[quad]);

    if (Mesh_conformize_cell_face(M, own, quad, M->nf, 4) != 0) return 1;

    E_Int nei = M->neigh[quad];

    if (nei != -1) {
        assert(M->clevel[nei] == M->flevel[quad]);

        if (Mesh_conformize_cell_face(M, nei, quad, M->nf, 4) != 0) return 1;
    }

    for (E_Int i = 0; i < 3; i++) {
        M->owner[M->nf+i] = own;
        M->neigh[M->nf+i] = nei;
    }

    // Update adaptation info
    M->flevel[quad]++;

    assert(M->fref[quad] == FACE_REFINED);

    for (E_Int i = 0; i < 3; i++) {
        E_Int fid = M->nf + i;

        M->flevel[fid] = M->flevel[quad];
        M->ftype[fid] = M->ftype[quad];

        M->fref[fid] = FACE_NEW;
    }

    M->fchildren[quad] = {quad, M->nf, M->nf+1, M->nf+2};

    M->fparent[M->nf]   = quad;
    M->fparent[M->nf+1] = quad;
    M->fparent[M->nf+2] = quad;

    M->ftag[M->nf]   = M->ftag[quad];
    M->ftag[M->nf+1] = M->ftag[quad];
    M->ftag[M->nf+2] = M->ftag[quad];

    //assert(M->fpattern[quad] == DIR_ISO);
    M->fpattern[M->nf] = M->fpattern[quad];
    M->fpattern[M->nf+1] = M->fpattern[quad];
    M->fpattern[M->nf+2] = M->fpattern[quad];

    // Increment face/edge/point count
    M->nf += 3;
    M->np += 1;

    return 0;
}

void Q4_reorder(E_Int *pn, E_Int reorient, E_Int i0, E_Int local[4])
{
    for (E_Int i = 0; i < 4; i++) local[i] = pn[i];
    Right_shift(local, i0, 4);
    if (reorient) std::swap(local[1], local[3]);
}

void refine_face_iso(E_Int face, Mesh *M)
{
    switch (M->ftype[face]) {
        case QUAD:
            Q9_refine(face, M);
            break;
        default:
            assert(0);
            break;
    }
}