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
#include "Hexa.h"
#include "Tetra.h"
#include "Penta.h"
#include "Pyra.h"
#include "Quad.h"
#include "Tri.h"

inline
E_Int Mesh_is_face_aligned_with_vec3(Mesh *M, E_Int fid, E_Float *vec3)
{
    E_Int pn[4];
    E_Int *face = Mesh_get_face(M, fid);
    E_Int *frange = Mesh_get_frange(M, fid);
    for (E_Int i = 0; i < M->fstride[fid]; i++) {
        pn[i] = face[2*i];
    }

    // Compute normal to the triangle formed by the first 3 points
    E_Float e0[3] = { M->X[pn[1]] - M->X[pn[0]],
                    M->Y[pn[1]] - M->Y[pn[0]],
                    M->Z[pn[1]] - M->Z[pn[0]] };

    E_Float e1[3] = { M->X[pn[2]] - M->X[pn[0]],
                    M->Y[pn[2]] - M->Y[pn[0]],
                    M->Z[pn[2]] - M->Z[pn[0]] };
    
    E_Float n[3];
    K_MATH::cross(e0, e1, n);

    E_Float dp = fabs(K_MATH::dot(n, vec3, 3)) / K_MATH::norm(n, 3);

    return fabs(dp - 1.0) <= 1e-8;
}

static inline
E_Int set_PENTA_for_2D(E_Int cid, Mesh *M)
{
    assert(0);
    return 1;
}

static inline
E_Int set_HEXA_for_2D(E_Int cid, Mesh *M)
{
    E_Int *cell = Mesh_get_cell(M, cid);
    E_Int *crange = Mesh_get_crange(M, cid);
    
    E_Int start = 0;
    for (; start < M->cstride[cid]; start++) {
        E_Int fid = cell[4*start];
        if (Mesh_is_face_aligned_with_vec3(M, fid, M->mode_2D)) break;
    }

    if (start == M->cstride[cid]) {
        merr("%d -> Cell " SF_D_ " is not aligned with 2D vector.", M->pid, cid);
        assert(0);
        return 1;
    }

    Right_shift(cell, 4*start, 24);
    Right_shift(crange, start, 6);

    return 0;
}

E_Int Mesh_set_cells_for_2D(Mesh *M)
{
    E_Int ret = 0;

    for (E_Int cid = 0; cid < M->nc; cid++) {
        E_Int ctype = M->ctype[cid];
        
        switch (ctype) {
            case HEXA:
                ret |= set_HEXA_for_2D(cid, M);
                break;
            case PENTA:
                assert(0);
                ret |= set_PENTA_for_2D(cid, M);
                break;
            default:
                merr("Mesh is not 2D.");
                assert(0);
                return 1;
        }
    }

    assert(ret == 0);
    return ret;
}

E_Int Mesh_set_face_types(Mesh *M)
{
    assert(M->ftype == NULL);

    M->ftype = IntArray(M->nf);

    for (E_Int i = 0; i < M->nf; i++) {
        E_Int *frange = Mesh_get_frange(M, i);
        E_Int np = 0;
        for (E_Int j = 0; j < M->fstride[i]; j++) np += frange[j];

        switch (np) {
            case 4:
                M->ftype[i] = QUAD;
                break;
            case 3:
                M->ftype[i] = TRI;
                break;
            default:
                fprintf(stderr, "Face " SF_D_ " is neither a TRI nor a QUAD.", i);
                assert(0);
                return 1;
        }
    }

    return 0;
}

E_Int Mesh_set_cell_types(Mesh *M)
{
    assert(M->ctype == NULL);

    M->ctype = IntArray(M->nc);

    for (E_Int i = 0; i < M->nc; i++) {
        E_Int nf = M->cstride[i];
        assert(nf == 6);

        switch (nf) {
            case 6:
                M->ctype[i] = HEXA;
                break;
            case 4:
                M->ctype[i] = TETRA;
                break;
            case 5: {
                // PENTA: 2 TRI + 3 QUAD
                // PYRA: 4 TRI + 1 QUAD
                E_Int ntri = 0, nquad = 0;

                E_Int *cell = Mesh_get_cell(M, i);

                for (E_Int j = 0; j < nf; j++) {
                    E_Int face = cell[4*j];
                    M->ftype[face] == TRI ? ntri++ : nquad++;
                }
                if (ntri == 2 && nquad == 3) {
                    M->ctype[i] = PENTA;
                } else if (ntri == 4 && nquad == 1) {
                    M->ctype[i] = PYRA;
                } else {
                    fprintf(stderr, "Cell " SF_D_ " has five faces but is "
                            "neither a PENTA nor a PYRA.\n", i);
                    assert(0);
                    return 1;
                }
                break;
            }
            default: {
                fprintf(stderr, "Cell " SF_D_ " is not one of the following"
                        "types: TETRA, HEXA, PENTA, PYRA.\n", i);
                assert(0);
                return 1;
            }
        }
    }
    
    return 0;
}
