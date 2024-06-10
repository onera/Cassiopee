#pragma once

#include "vec3.h"

typedef Vec3 point;

struct pointFace {
    E_Int F;
    E_Int T;

    pointFace(E_Int face, E_Int tri);
};