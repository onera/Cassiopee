#pragma once

#include <vector>

#include "point.h"
#include "smesh.h"

struct Triangle {
    E_Int a, b, c;

    static E_Int ispointInside(E_Float x, E_Float y, E_Float ax, E_Float ay,
        E_Float bx, E_Float by, E_Float cx, E_Float cy);

    static E_Int ray_E_Intersect(E_Float px, E_Float py, E_Float pz, E_Float dx,
        E_Float dy, E_Float dz, E_Float ax, E_Float ay, E_Float az, E_Float bx,
        E_Float by, E_Float bz, E_Float cx, E_Float cy, E_Float cz, E_Float &u,
        E_Float &v, E_Float &w, E_Float &t, E_Float &x, E_Float &y, E_Float &z);
};

struct TriMesh {
    std::vector<point> P;
    std::vector<Triangle> T;

    TriMesh(const Smesh &M);

    pointFace locate(point p);
};
