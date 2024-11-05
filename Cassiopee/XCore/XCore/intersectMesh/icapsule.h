#pragma once

#include "mesh.h"

struct ICapsule {
    IMesh M;
    std::vector<IMesh> Ss;

    ICapsule(const Karray &marray, const std::vector<Karray> &sarrays,
        const std::vector<E_Float *> &ptags);
    
    static std::vector<PointLoc> refine(Smesh &Mf,
        std::set<E_Int> &mfids, Smesh &Sf);
};
