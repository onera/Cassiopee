/*    
    Copyright 2013-2026 ONERA.

    This file is part of Cassiopee.

    Cassiopee is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    Cassiopee is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Cassiopee.  If not, see <http://www.gnu.org/licenses/>.
*/

#pragma once

#include "xcore.h"
#include <vector>

struct CompensatedSum {
    E_Float x, y;

    CompensatedSum &operator-()
    {
        x = -x;
        y = -y;
        return *this;
    }
};

struct Expansion {
    std::vector<E_Float> m_terms;

    explicit Expansion(const CompensatedSum &s)
    {
        m_terms.resize(2);
        m_terms[0] = s.y;
        m_terms[1] = s.x;
    }

    explicit Expansion(size_t m) { m_terms.resize(m); }

    size_t size() const { return m_terms.size(); }

    E_Float operator[](E_Int idx) const
    {
        assert(idx >= 0);
        assert((size_t)idx < m_terms.size());
        return m_terms[idx];
    }

    E_Float &operator[](E_Int idx)
    {
        assert(idx >= 0);
        assert((size_t)idx < m_terms.size());
        return m_terms[idx];
    }

    void print() const
    {
        for (auto term : m_terms) {
            printf("%g + ", term);
        }
        printf("\n");
    }

    int sign() const;
};

CompensatedSum two_sum(E_Float a, E_Float b);
CompensatedSum fast_two_sum(E_Float a, E_Float b);
CompensatedSum two_diff(E_Float a, E_Float b);
CompensatedSum two_product(E_Float a, E_Float b);
Expansion linear_sum(const Expansion &e, const Expansion &f);
CompensatedSum split(E_Float a);
Expansion difference_of_products(E_Float a, E_Float b, E_Float c, E_Float d);