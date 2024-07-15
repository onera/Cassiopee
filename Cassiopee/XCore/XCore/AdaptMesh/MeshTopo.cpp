#include "Mesh.h"
#include "../common/mem.h"
#include "Hexa.h"
#include "Tetra.h"
#include "Penta.h"
#include "Pyra.h"
#include "Quad.h"
#include "Tri.h"

inline
Int Mesh_is_face_aligned_with_vec3(Mesh *M, Int fid, Float *vec3)
{
    Int nn, pn[8];
    Mesh_get_fpoints(M, fid, nn, pn);

    // Compute normal to the triangle formed by the first 3 points
    Float e0[3] = { M->X[pn[1]] - M->X[pn[0]],
                    M->Y[pn[1]] - M->Y[pn[0]],
                    M->Z[pn[1]] - M->Z[pn[0]] };

    Float e1[3] = { M->X[pn[2]] - M->X[pn[0]],
                    M->Y[pn[2]] - M->Y[pn[0]],
                    M->Z[pn[2]] - M->Z[pn[0]] };
    
    Float n[3];
    K_MATH::cross(e0, e1, n);

    Float dp = fabs(K_MATH::dot(n, vec3, 3)) / K_MATH::norm(n, 3);

    return fabs(dp - 1.0) <= 1e-8;
}

static inline
Int set_PENTA_for_2D(Int cid, Mesh *M)
{
    assert(0);
    return 1;
}

static inline
Int set_HEXA_for_2D(Int cid, Mesh *M)
{
    Int *cell = Mesh_get_cell(M, cid);
    Int *crange = Mesh_get_crange(M, cid);
    
    Int start = 0;
    for (; start < M->cstride[cid]; start++) {
        Int fid = cell[4*start];
        if (Mesh_is_face_aligned_with_vec3(M, fid, M->mode_2D)) break;
    }

    if (start == M->cstride[cid]) {
        merr("Cell %d is not aligned with 2D vector.", cid);
        assert(0);
        return 1;
    }

    Right_shift(cell, 4*start, 24);
    Right_shift(crange, start, 6);

    return 0;
}

Int Mesh_set_cells_for_2D(Mesh *M)
{
    Int ret = 0;

    for (Int cid = 0; cid < M->nc; cid++) {
        Int ctype = M->ctype[cid];
        
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

Int Mesh_set_face_types(Mesh *M)
{
    assert(M->ftype == NULL);

    M->ftype = IntArray(M->nf);

    for (Int i = 0; i < M->nf; i++) {
        Int *frange = Mesh_get_frange(M, i);
        Int np = 0;
        for (Int j = 0; j < M->fstride[i]; j++) np += frange[j];

        switch (np) {
            case 4:
                M->ftype[i] = QUAD;
                break;
            case 3:
                M->ftype[i] = TRI;
                break;
            default:
                fprintf(stderr, "Face %d is neither a TRI nor a QUAD.", i);
                assert(0);
                return 1;
        }
    }

    return 0;
}

Int Mesh_set_cell_types(Mesh *M)
{
    assert(M->ctype == NULL);

    M->ctype = IntArray(M->nc);

    for (Int i = 0; i < M->nc; i++) {
        Int nf = M->cstride[i];
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
                Int ntri = 0, nquad = 0;

                Int *cell = Mesh_get_cell(M, i);

                for (Int j = 0; j < nf; j++) {
                    Int face = cell[4*j];
                    M->ftype[face] == TRI ? ntri++ : nquad++;
                }
                if (ntri == 2 && nquad == 3) {
                    M->ctype[i] = PENTA;
                } else if (ntri == 4 && nquad == 1) {
                    M->ctype[i] = PYRA;
                } else {
                    fprintf(stderr, "Cell %d has five faces but is neither a"
                            " PENTA nor a PYRA.\n", i);
                    assert(0);
                    return 1;
                }
                break;
            }
            default: {
                fprintf(stderr, "Cell %d is not one of the following types"
                            ": TETRA, HEXA, PENTA, PYRA.\n", i);
                assert(0);
                return 1;
            }
        }
    }
    
    return 0;
}
