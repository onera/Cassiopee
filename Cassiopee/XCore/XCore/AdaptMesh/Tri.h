#pragma once

#include "common/common.h"
#include "Mesh.h"

void refine_tri(Int tri, Mesh *M);

inline void T6_get_ordered_data(Mesh *M, Int node, Int reorient, Int *children,
    Int local[6])
{
}

void reorder_tri(Int tri, Mesh *M);

Int check_canon_tri(Int tri, Mesh *M);