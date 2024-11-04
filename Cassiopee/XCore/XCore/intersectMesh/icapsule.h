#pragma once

#include "mesh.h"

struct ICapsule {
    IMesh M;
    std::vector<IMesh> Ss;

    static ICapsule do_it(const Karray &marray,
        const std::vector<Karray> &sarray,
        const std::vector<E_Float *> &ptags);
    
    static std::vector<PointLoc> refine(Smesh &Mf,
        std::set<E_Int> &mfids, Smesh &Sf);
};
