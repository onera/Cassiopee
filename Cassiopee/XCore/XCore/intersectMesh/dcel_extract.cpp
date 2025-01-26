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
#include "dcel.h"
#include "smesh.h"

Smesh Dcel::export_smesh(bool check_Euler) const
{
    Smesh smesh;
    smesh.check_Euler = check_Euler;

    smesh.np = V.size();
    smesh.X.resize(smesh.np);
    smesh.Y.resize(smesh.np);
    smesh.Z.resize(smesh.np);

    for (Vertex *v : V) {
        smesh.X[v->id] = v->x;
        smesh.Y[v->id] = v->y;
        smesh.Z[v->id] = v->z;
    }

    for (size_t i = 0; i < F.size(); i++) {
        Face *f = F[i];
        Hedge *REP = f->rep;
        Cycle *c = REP->cycle;
        if (c->inout != Cycle::OUTER) continue;
        Hedge *h = REP;

        std::vector<E_Int> pn;

        do {
            Vertex *v = h->orig;
            pn.push_back(v->id);
            h = h->next;
        } while (h != REP);

        smesh.F.push_back(pn);
    }

    smesh.nf = (E_Int)smesh.F.size();

    smesh.Fc = smesh.F;

    smesh.make_edges();
    
    return smesh;
}

std::vector<Dcel::Cycle *> Dcel::extract_cycles_of_indices(
    const std::vector<E_Int> &indices) const
{
    std::vector<Cycle *> ret;
    ret.reserve(indices.size());

    for (E_Int index : indices) {
        ret.push_back(C[index]);
    }

    return ret;
}

std::vector<E_Int> Dcel::extract_indices_of_type(int type) const
{
    std::vector<E_Int> ret;

    for (size_t i = 0; i < C.size(); i++) {
        if (C[i]->inout == type) ret.push_back(i);
    }

    return ret;
}
