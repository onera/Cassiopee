/*    
    Copyright 2013-2024 Onera.

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
#include <cstdio>

#include "smesh.h"

void Smesh::refine_edge(const u_edge &e)
{
    E_Int p = e.p;
    E_Int q = e.q;

    X[np] = 0.5 * (X[p] + X[q]);
    Y[np] = 0.5 * (Y[p] + Y[q]);
    Z[np] = 0.5 * (Z[p] + Z[q]);

    ecenter[e] = np;
    np++;
}

void Smesh::refine_tri(E_Int fid)
{
    std::vector<E_Int> nodes(F[fid]);

    E_Int ec[3];

    for (size_t i = 0; i < nodes.size(); i++) {
        E_Int p = nodes[i];
        E_Int q = nodes[(i+1)%nodes.size()];
        u_edge e(p, q);
        auto it = ecenter.find(e);
        if (it == ecenter.end()) {
            refine_edge(e);
            ec[i] = ecenter[e];
        } else {
            ec[i] = it->second;
        }
    }

    fchildren[fid].push_back({nf, nf+1, nf+2});

    F[fid]  = {nodes[0], ec[0], ec[2]};
    F[nf]   = {ec[0], nodes[1], ec[1]};
    F[nf+1] = {ec[2], ec[1], nodes[2]};
    F[nf+2] = {ec[0], ec[1], ec[2]};
}

void Smesh::resize_for_refinement(size_t nref_faces)
{
    E_Int fincr = nref_faces * 3;
    E_Int pincr = nref_faces * 3;
    F.resize(nf + fincr);
    X.resize(np + pincr);
    Y.resize(np + pincr);
    Z.resize(np + pincr);
}

void Smesh::get_leaves(E_Int fid, std::vector<E_Int> &leaves) const
{
}
