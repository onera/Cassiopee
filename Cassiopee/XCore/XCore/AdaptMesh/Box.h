#pragma once

#include "xcore.h"

struct Mesh;
struct FaceSort;
struct DynMesh;

struct Box3 {
    E_Float xmin, ymin, zmin;
    E_Float xmax, ymax, zmax;
};

inline
bool Box3_in_Box3(const Box3 small, const Box3 big)
{
    return (big.xmin <= small.xmin) && (big.xmax >= small.xmax) &&
           (big.ymin <= small.ymin) && (big.ymax >= small.ymax) &&
           (big.zmin <= small.zmin) && (big.zmax >= small.zmax);
}

Box3 Box3_make
(
    const Mesh *M,
    const FaceSort *mfaces,
    E_Int start, E_Int end
);

Box3 Box3_make
(
    const Mesh *M,
    const E_Int *skin,
    const E_Int *indices,
    E_Int start, E_Int end
);

void Box3_clamp(const Box3 *parent, Box3 *child);

/* DynMesh */

Box3 Box3_make
(
    const DynMesh *M,
    const FaceSort *mfaces,
    E_Int start, E_Int end
);

Box3 Box3_make
(
    const DynMesh *M,
    const E_Int *skin,
    const E_Int *indices,
    E_Int start, E_Int end
);