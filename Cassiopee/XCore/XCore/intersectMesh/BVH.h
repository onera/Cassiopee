#pragma once

#include "AABB.h"

struct BVH_node {
    AABB box = AABB_HUGE;
    E_Int left_node, first_tri_idx, tri_count;
    bool is_leaf() const { return tri_count > 0; }
};
