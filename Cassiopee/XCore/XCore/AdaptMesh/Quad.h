#pragma once

#include "common/common.h"
#include "Mesh.h"

void Q9_refine(Int quad, Mesh *M);

void Q6_refine(Int quad, Mesh *M);

void Q4_reorder(Int *pn, Int reorient, Int i0, Int local[4]);

inline
void refine_face_iso(Int face, Mesh *M)
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