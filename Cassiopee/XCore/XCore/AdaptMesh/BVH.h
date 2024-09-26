#pragma once

#include "Box.h"

struct Mesh;
struct FaceSort;
struct Vec3f;
struct Point;
struct PointFaces;
struct DynMesh;

struct BVH_node {
    Box3 box;
    E_Int start, end;
    BVH_node *left;
    BVH_node *right;
};

BVH_node *BVH_create_node(const Box3 *box, E_Int start, E_Int end,
    BVH_node *left, BVH_node *right);

BVH_node *BVH_make(const Mesh *M, FaceSort *mfaces, E_Int start, E_Int end,
    const Box3 *parent_box);

BVH_node *BVH_make(const Mesh *M, const E_Int *skin, const Vec3f *fc,
    E_Int *indices, E_Int start, E_Int end, const Box3 *parent_box);

void BVH_locate_point
(
    const BVH_node *node,
    const Mesh *M,
    const E_Int *skin,
    const E_Int *indices,
    const Point *p,
    PointFaces *pfaces
);

void BVH_free(BVH_node *node);

/* DynMesh */

BVH_node *BVH_make(const DynMesh *M, FaceSort *mfaces, E_Int start, E_Int end,
    const Box3 *parent_box);

BVH_node *BVH_make(const DynMesh *M, const E_Int *skin, const Vec3f *fc,
    E_Int *indices, E_Int start, E_Int end, const Box3 *parent_box);

void BVH_locate_point
(
    const BVH_node *node,
    const DynMesh *M,
    const E_Int *skin,
    const E_Int *indices,
    const Point *p,
    PointFaces *pfaces
);