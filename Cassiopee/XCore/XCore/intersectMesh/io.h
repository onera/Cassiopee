#pragma once

#include <vector>

#include "xcore.h"
#include "vertex.h"

void point_write(const char *fname, const std::vector<Vertex *> &I);

void point_write(const char *fname, E_Float *Xs, E_Float *Ys, E_Float *Zs,
    const std::vector<E_Int> &proj_points);