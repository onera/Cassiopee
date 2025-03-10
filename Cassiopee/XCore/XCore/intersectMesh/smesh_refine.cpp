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
#include "smesh.h"

void Smesh::update_plocs(const std::vector<E_Int> &fparents,
    std::vector<PointLoc> &plocs)
{
    // Children to update: fchildren[fparent].back();
    assert(0);
}

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

inline
void Smesh::compute_face_center(E_Int fid)
{
    E_Float *fc = &fcenters[3*fid];
    fc[0] = fc[1] = fc[2] = 0;
    const auto &pn = Fc[fid];
    for (E_Int p : pn) {
        fc[0] += X[p];
        fc[1] += Y[p];
        fc[2] += Z[p];
    }
    for (E_Int i = 0; i < 3; i++) fc[i] /= pn.size();
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

    Fc[fid]  = F[fid];
    Fc[nf]   = F[nf];
    Fc[nf+1] = F[nf+1];
    Fc[nf+2] = F[nf+2];

    const E_Float *N = &fnormals[3*fid];
    for (E_Int i = 0; i < 3; i++) {
        E_Float *fN = &fnormals[3*(nf+i)];
        for (E_Int j = 0; j < 3; j++) fN[j] = N[j];
    }

    //compute_face_center(fid);
    //for (E_Int i = 0; i < 3; i++) compute_face_center(nf+i);

    nf += 3;
}

void Smesh::resize_for_refinement(size_t nref_faces)
{
    E_Int fincr = nref_faces * 3;
    E_Int pincr = nref_faces * 3;
    F.resize(nf + fincr);
    Fc.resize(nf + fincr);
    X.resize(np + pincr, EFLOATMAX);
    Y.resize(np + pincr, EFLOATMAX);
    Z.resize(np + pincr, EFLOATMAX);
    fnormals.resize(3*(nf + fincr), EFLOATMAX);
    fcenters.resize(3*(nf + fincr), EFLOATMAX);
}

void Smesh::refine(const std::vector<E_Int> &ref_faces)
{
    resize_for_refinement(ref_faces.size());
    for (E_Int fid : ref_faces)
        refine_tri(fid);
}
