#pragma once

#include "common/common.h"
#include "Mesh.h"

const Int normalIn_H[6] = {1, 0, 1, 0, 1, 0};

void H27_refine(Int hexa, Mesh *M);

void H27_reorder(Int hexa, Mesh *M);

void H18_refine(Int hexa, Mesh *M);

void H18_reorder(Int hexa, Mesh *M);

Int check_canon_hexa(Int hexa, Mesh *M);

inline
void update_shell_pe(Int hexa, Mesh *M)
{
    const auto &children = M->cchildren.at(hexa);

    for (Int cid : children) {
        Int *child = Mesh_get_cell(M, cid);

        for (Int j = 0; j < 6; j++) {
            Int face = child[4*j];
            
            if      (M->owner[face] == hexa) M->owner[face] = cid;
            else if (M->neigh[face] == hexa) M->neigh[face] = cid;
        }
    }
}

inline
void update_range_and_stride(Mesh *M, Int hexa, Int cpos, Int nchildren)
{
    Int *crange = Mesh_get_crange(M, hexa);
    for (Int i = 0; i < M->cstride[hexa]; i++) {
        crange[i] = 1;
    }

    for (Int i = 0; i < nchildren; i++) {
        Int child = cpos + i;

        M->cstride[child] = M->cstride[hexa];

        crange = Mesh_get_crange(M, child);
        for (Int j = 0; j < M->cstride[child]; j++) {
            crange[j] = 1;
        }
    }
}