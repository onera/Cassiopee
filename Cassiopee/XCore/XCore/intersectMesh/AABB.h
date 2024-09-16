#pragma once

#include "xcore.h"

struct IMesh;

struct AABB {
    E_Float xmin;
    E_Float ymin;
    E_Float zmin;
    E_Float xmax;
    E_Float ymax;
    E_Float zmax;

    AABB();

    AABB(const IMesh &M, E_Int *ids, E_Int count);
};