#pragma once

#include "xcore.h"

struct u_edge {
    E_Int p, q;

    u_edge(E_Int P, E_Int Q)
    : p(std::min(P, Q)), q(std::max(P, Q))
    {}

    bool operator<(const u_edge& e) const
    {
        return (p < e.p) || (p == e.p && q < e.q);
    }
};