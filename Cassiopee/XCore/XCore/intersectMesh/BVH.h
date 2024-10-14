#pragma once

#include "AABB.h"

struct BVH_node {
    AABB box = AABB_HUGE;
    E_Int start = -1;
    E_Int end = -1;
    BVH_node *left = NULL;
    BVH_node *right = NULL;
};