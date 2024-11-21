#pragma once

#include "xcore.h"

struct Mesh;

struct FaceSort {
    E_Int fid;
    E_Float fc[3];
    E_Float UX, UY, UZ;
    E_Float VX, VY, VZ;
    E_Float UU, VV, UV;
    E_Float inv_denom;
    E_Float xa, ya, za;
};

void FaceSort_compute_data(const Mesh *M, FaceSort *mfaces, E_Int mcount);