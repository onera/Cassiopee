/*    
    Copyright 2013-2025 Onera.

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

#include "precise.h"

static
std::vector<E_Float> merge_and_sort(const Expansion &e, const Expansion &f)
{
    std::vector<E_Float> sequence;
    sequence.reserve(e.size() + f.size());
    for (size_t i = 0; i < e.size(); i++) sequence.push_back(e[i]);
    for (size_t i = 0; i < f.size(); i++) sequence.push_back(f[i]);
    std::sort(sequence.begin(), sequence.end(),
        [&] (E_Float a, E_Float b)
        {
            return fabs(a) < fabs(b);
        }
    );
    return sequence;
}

CompensatedSum fast_two_sum(E_Float a, E_Float b)
{
    assert(fabs(a) >= fabs(b));
    E_Float x = a + b;
    E_Float bvirtual = x - a;
    E_Float y = b - bvirtual;
    return {x, y};
}

CompensatedSum two_sum(E_Float a, E_Float b)
{
    E_Float x = a + b;
    E_Float bvirtual = x - a;
    E_Float avirtual = x - bvirtual;
    E_Float broundoff = b - bvirtual;
    E_Float aroundoff = a - avirtual;
    E_Float y = aroundoff + broundoff;
    return {x, y};
}

CompensatedSum two_diff(E_Float a, E_Float b)
{
    E_Float x = a - b;
    E_Float bvirtual = a - x;
    E_Float avirtual = x + bvirtual;
    E_Float broundoff = bvirtual - b;
    E_Float aroundoff = a - avirtual;
    E_Float y = aroundoff + broundoff;
    return {x, y};
}

Expansion linear_sum(const Expansion &e, const Expansion &f)
{
    size_t m = e.size(), n = f.size();
    Expansion h(m + n);
    auto g = merge_and_sort(e, f);
    CompensatedSum Q = fast_two_sum(g[1], g[0]);
    for (int i = 2; i < m + n; i++) {
        CompensatedSum R = fast_two_sum(g[i], Q.y);
        h[i-2] = R.y;
        Q = two_sum(Q.x, R.x);
    }
    h[m+n-2] = Q.y;
    h[m+n-1] = Q.x;
    return h;
}

// Note(Imad): for the 64bit double of C
constexpr E_Float splitter = (1 << 27) + 1;

CompensatedSum split(E_Float a)
{
    E_Float c = splitter * a;
    E_Float abig = c - a;
    E_Float ahi = c - abig;
    E_Float alo = a - ahi;
    return {ahi, alo};
}

CompensatedSum two_product(E_Float a, E_Float b)
{
    E_Float x = a * b;
    CompensatedSum A = split(a);
    CompensatedSum B = split(b);
    E_Float err1 = x - (A.x * B.x);
    E_Float err2 = err1 - (A.y * B.x);
    E_Float err3 = err2 - (A.x * B.y);
    E_Float y = (A.y * B.y) - err3;
    return {x, y};
}

Expansion difference_of_products(E_Float a, E_Float b, E_Float c, E_Float d)
{
    auto S1 = two_product(a, b);
    auto S2 = -two_product(c, d);
    Expansion E1(S1);
    Expansion E2(S2);
    Expansion E = linear_sum(E1, E2);
    return E;
}

int Expansion::sign() const
{
    for (E_Int i = m_terms.size()-1; i >= 0; i--) {
        if (m_terms[i] == 0.0) continue;
        if (m_terms[i] < 0.0) return -1;
        return 1;
    }
    return 0;
}
