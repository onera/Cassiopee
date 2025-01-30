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
#include "Quad.h"
#include "Mesh.h"

E_Int Q6_refine(E_Int quad, Mesh *M)
{
    E_Int NODES[6];
    for (E_Int i = 0; i < 6; i++) NODES[i] = -1;

    E_Int *fpts = Mesh_get_face(M, quad);
    E_Int *frange = Mesh_get_frange(M, quad);

    NODES[0] = fpts[0];
    NODES[1] = fpts[2];
    NODES[2] = fpts[4];
    NODES[3] = fpts[6];

    // Get the first edge that's perpendicular to 2D normal vector

    E_Int start = -1;
    for (E_Int i = 0; i < 8; i += 2) {
        E_Int p = fpts[i];
        E_Int q = fpts[(i+2)%8];

        E_Float e[3] = { M->X[q] - M->X[p],
                       M->Y[q] - M->Y[p],
                       M->Z[q] - M->Z[p] };
        
        E_Float dp = K_MATH::dot(e, M->mode_2D, 3);

        if (fabs(dp) <= 1e-8) {
            start = i;
            break;
        }
    }

    assert(start != -1);

    E_Int idx[8];
    for (E_Int i = 0; i < 8; i++) idx[i] = i;

    E_Int i0 = Get_pos(start, idx, 8);
    Right_shift(idx, i0, 8);
    assert(idx[0] == start);

    // Refine the first edge
    if (frange[idx[0]/2] == 2) {
        NODES[4] = fpts[idx[1]];
    } else {
        E_Int p = fpts[idx[0]];
        E_Int q = fpts[idx[2]];
        Mesh_refine_or_get_edge_center(M, p, q, NODES[4]);
    }

    // Refine the second edge
    if (frange[idx[4]/2] == 2) {
        NODES[5] = fpts[idx[5]];
    } else {
        E_Int p = fpts[idx[4]];
        E_Int q = fpts[idx[6]];
        Mesh_refine_or_get_edge_center(M, p, q, NODES[5]);
    }

    // Setup face points

    if (start == 0 || start == 4) {
        // Dir X

        fpts = Mesh_get_face(M, M->nf);
        fpts[0] = NODES[4]; fpts[2] = NODES[1];
        fpts[4] = NODES[2]; fpts[6] = NODES[5];

        fpts = Mesh_get_face(M, quad);
        memset(fpts, -1, 8*sizeof(E_Int));
        fpts[0] = NODES[0]; fpts[2] = NODES[4];
        fpts[4] = NODES[5]; fpts[6] = NODES[3];
    } else {
        assert(start == 2 || start == 6);

        // Dir Y

        fpts = Mesh_get_face(M, M->nf);
        fpts[0] = NODES[5]; fpts[2] = NODES[4];
        fpts[4] = NODES[2]; fpts[6] = NODES[3];

        fpts = Mesh_get_face(M, quad);
        memset(fpts, -1, 8*sizeof(E_Int));
        fpts[0] = NODES[0]; fpts[2] = NODES[1];
        fpts[4] = NODES[4]; fpts[6] = NODES[5];
    }

    // Update ranges and strides
    Mesh_update_face_range_and_stride(M, quad, M->nf, 1);

    // Conformize parent cells

    E_Int own = M->owner[quad];

    assert(M->clevel[own] == M->flevel[quad]);

    if (Mesh_conformize_cell_face(M, own, quad, M->nf, 2) != 0) return 1;

    E_Int nei = M->neigh[quad];

    if (nei != -1) {
        assert(M->clevel[nei] == M->flevel[quad]);

        if (Mesh_conformize_cell_face(M, nei, quad, M->nf, 2) != 0) return 1;
    }

    for (E_Int i = 0; i < 1; i++) {
        M->owner[M->nf+i] = own;
        M->neigh[M->nf+i] = nei;
    }

    // Update adaptation info
    M->flevel[quad]++;

    assert(M->fref[quad] == FACE_REFINED);

    for (E_Int i = 0; i < 1; i++) {
        E_Int fid = M->nf + i;

        M->flevel[fid] = M->flevel[quad];
        M->ftype[fid] = M->ftype[quad];

        M->fref[fid] = FACE_NEW;
    }

    M->fchildren[quad] = {quad, M->nf};

    M->fparent[M->nf] = quad;

    M->ftag[M->nf] = M->ftag[quad];

    // Increment face/edge/point count
    M->nf += 1;

    return 0;
}
