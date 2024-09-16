#pragma once

#include "AABB.h"
#include "common/mem.h"

#define MAX_FACES_PER_LEAF 10

struct IMesh;

struct BVH_node {
    AABB bbox;
    BVH_node *left;
    BVH_node *right;
    E_Int *elements;
};

BVH_node *BVH_create_node(AABB aabb, BVH_node *left_child,
    BVH_node *right_child, Mem_arena &arena);

BVH_node *BVH_create_node(AABB aabb, BVH_node *left_child, BVH_node *right_child,
        E_Int *ids, E_Int count, Mem_arena &arena);

BVH_node *BVH_create_node(E_Int *ids, E_Int count, const IMesh &M,
    Mem_arena &arena);

struct BVH {
    BVH_node *root;

    Mem_arena arena;

    BVH(const IMesh &M);
};