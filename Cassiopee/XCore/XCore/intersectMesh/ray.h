#pragma once

#include "triangleIntersection.h"

E_Int MollerTrumbore(E_Float px, E_Float py, E_Float pz, E_Float dx, E_Float dy,
    E_Float dz, E_Float ax, E_Float ay, E_Float az, E_Float bx, E_Float by, E_Float bz,
    E_Float cx, E_Float cy, E_Float cz, TriangleIntersection &TI);