#pragma once

#include <vector>
#include <unordered_map>

#include "xcore.h"
#include "vertex.h"
#include "triangleIntersection.h"

void point_write(const char *fname, const std::vector<Vertex *> &I);

void point_write(const char *fname, E_Float *Xs, E_Float *Ys, E_Float *Zs,
    const std::vector<E_Int> &proj_points);

void edge_write(const char *fname, E_Float *X, E_Float *Y, E_Float *Z,
    const std::unordered_map<E_Int, TriangleIntersection> &point_hits);