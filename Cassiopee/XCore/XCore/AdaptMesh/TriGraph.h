#pragma once

#include "common/mem.h"
#include <vector>
#include <array>

struct TriGraph {
    size_t nf;
    std::vector<E_Int> skin;
    std::vector<std::array<E_Int, 3>> T;
    std::vector<std::array<E_Int, 3>> E;
    std::vector<E_Int> fdat;
    std::vector<E_Int> level;
    
    void clear();
}; 


