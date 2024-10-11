#pragma once

#include "mesh.h"

struct ICapsule {
    IMesh M;
    std::vector<IMesh> Ss;
    std::vector<Smesh> mpatches;
    std::vector<Smesh> spatches;

    ICapsule(const Karray &marray, const std::vector<Karray> &sarray,
        const std::vector<E_Float *> &ptags);
};
