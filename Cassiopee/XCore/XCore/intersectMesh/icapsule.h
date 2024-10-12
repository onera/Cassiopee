#pragma once

#include "mesh.h"

struct ICapsule {
    IMesh M;
    std::vector<IMesh> Ss;
    std::vector<Smesh> mpatches;
    std::vector<Smesh> spatches;

    E_Float NEAR_VERTEX_TOL = 1e-3;
    E_Float NEAR_EDGE_TOL = 1e-3;

    ICapsule(const Karray &marray, const std::vector<Karray> &sarray,
        const std::vector<E_Float *> &ptags);

    void correct_near_points_and_edges(Smesh &Sf, std::vector<PointLoc> &plocs,
        const IMesh &M);
};
