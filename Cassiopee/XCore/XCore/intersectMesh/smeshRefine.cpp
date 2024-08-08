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

// Refine wrt M's point cloud
size_t Smesh::refine(Smesh &M)
{
    // Locate M points within my faces
    std::vector<std::vector<pointFace>> pfs;
    pfs.reserve(M.np);

    for (Int i = 0; i < M.np; i++) {
        auto pf = locate(M.X[i], M.Y[i], M.Z[i]);
        pfs.push_back(pf);
    }

    // Face to mpoints located within it
    std::map<Int, std::vector<Int>> fpoints;

    for (size_t i = 0; i < pfs.size(); i++) {
        for (const auto &faces : pfs[i]) {
            Int F = faces.F;
            fpoints[F].push_back(i);
        }
    }

    // Keep the faces which contain 3 points or more
    std::map<Int, std::vector<Int>> filtered_faces_map;
    std::vector<Int> filtered_faces;

    for (auto &fpts : fpoints) {
        if (fpts.second.size() >= 3) {
            filtered_faces_map.insert({fpts.first, fpts.second});
            filtered_faces.push_back(fpts.first);
        }
    }

    // Make point-to-face connectivity for M
    M.make_point_faces();

    std::map<Int, std::vector<Int>> ref_faces_to_Mfaces;

    for (const auto &fdata : filtered_faces_map) {
        Int face = fdata.first;
        const auto &points = fdata.second;

        // For every mface, how many points are contained within face?
        std::map<Int, size_t> mcontained;

        for (Int p : points) {
            for (Int mface : M.P2F[p]) {
                assert(M.face_is_active(mface));
                mcontained[mface]++;
            }
        }

        for (const auto &mfdat : mcontained) {
            Int mface = mfdat.first;
            size_t npts = mfdat.second;

            assert(npts <= M.F[mface].size());
            
            // Do not keep mface if it is not completely enclosed by face
            if (npts < M.F[mface].size()) continue;

            // Do not keep mface if it is a duplicate of face
            if (faces_are_dups(face, mface, M)) continue;

            ref_faces_to_Mfaces[face].push_back(mface);
        }
    }
    
    printf("Faces to refine: %zu\n", ref_faces_to_Mfaces.size());

    filtered_faces.clear();

    for (auto &fdata : ref_faces_to_Mfaces) {
        filtered_faces.push_back(fdata.first);
    }

    Int iter = 0;
    size_t ret = 0;

    while (!ref_faces_to_Mfaces.empty()) {
        iter++;

        auto ref_data = smooth_ref_data(ref_faces_to_Mfaces);

        auto ref_faces = prepare_for_refinement(ref_data);

        ret += ref_faces.size();

        refine_faces(ref_faces);

        std::map<Int, std::vector<Int>> new_ref_faces_to_Mfaces;

        for (const auto &fdata : ref_faces_to_Mfaces) {
            Int parent = fdata.first;
            const auto &Mfaces = fdata.second;

            const auto &children = fchildren[parent];

            for (auto child : children) {
                for (Int mface : Mfaces) {
                    if (face_contains_Mface(child, mface, M)) {
                        assert(!faces_are_dups(child, mface, M));
                        new_ref_faces_to_Mfaces[child].push_back(mface);
                    }
                }
            }
        }

        ref_faces_to_Mfaces = new_ref_faces_to_Mfaces;
    }

    return ret;
}

void Smesh::refine_edge(Int edge)
{
    Int p = E[edge].p;
    Int q = E[edge].q;

    // Rounding errors!
    X[np] = 0.5 * (X[p] + X[q]);
    Y[np] = 0.5 * (Y[p] + Y[q]);

    Int ne1 = ne + 1;

    E[ne].p = p;
    E[ne].q = np;
    E[ne1].p = np;
    E[ne1].q = q;

    E2F[ne][0] = E2F[ne1][0] = E2F[edge][0];
    E2F[ne][1] = E2F[ne1][1] = E2F[edge][1];

    eactive.erase(edge);
    eactive.insert(ne);
    eactive.insert(ne1);

    echildren[edge] = {ne, ne1};
    
    elevel[ne] = elevel[ne+1] = elevel[edge] + 1;

    l2gp[np] = M_np++;

    l2ge[ne] = M_ne++;
    l2ge[ne+1] = M_ne++;

    np += 1;
    ne += 2;
}

void Smesh::refine_quad(Int quad)
{
    // Refine the edges
    const auto &pe = F2E[quad];
    assert(pe.size() == 4);
    for (Int edge : pe) {
        if (echildren.find(edge) == echildren.end()) {
            refine_edge(edge);
        }
    }

    // Edge centers
    Int ec[4];
    for (Int i = 0; i < 4; i++) {
        ec[i] = get_edge_center(pe[i]);
    }

    // Face centroid
    const auto &pn = F[quad];

    X[np] = Y[np] = 0.0;
    for (Int i = 0; i < 4; i++) {
        Int p = pn[i];
        X[np] += X[p];
        Y[np] += Y[p];
    }
    X[np] *= 0.25;
    Y[np] *= 0.25;

    // New face points
    Int nf0 = nf, nf1 = nf+1, nf2 = nf+2, nf3 = nf+3;

    F[nf0] = { pn[0], ec[0], np,    ec[3] };
    F[nf1] = { ec[0], pn[1], ec[1], np    };
    F[nf2] = { np   , ec[1], pn[2], ec[2] };
    F[nf3] = { ec[3], np   , ec[2], pn[3] };

    // Internal edge points
    Int ne0 = ne, ne1 = ne+1, ne2 = ne+2, ne3 = ne+3;
    
    E[ne0].p = ec[0]; E[ne0].q = np;
    E[ne1].p = ec[1]; E[ne1].q = np;
    E[ne2].p = ec[2]; E[ne2].q = np;
    E[ne3].p = ec[3]; E[ne3].q = np;

    // New face edges
    const auto &pe0 = echildren[pe[0]];
    const auto &pe1 = echildren[pe[1]];
    const auto &pe2 = echildren[pe[2]];
    const auto &pe3 = echildren[pe[3]];

    Int eid[8];
    eid[0] = pe0[0]; eid[1] = pe0[1];
    eid[2] = pe1[0]; eid[3] = pe1[1];
    eid[4] = pe2[0]; eid[5] = pe2[1];
    eid[6] = pe3[0]; eid[7] = pe3[1];

    if (E2F[pe[0]][0] != quad) std::swap(eid[0], eid[1]);
    if (E2F[pe[1]][0] != quad) std::swap(eid[2], eid[3]);
    if (E2F[pe[2]][0] != quad) std::swap(eid[4], eid[5]);
    if (E2F[pe[3]][0] != quad) std::swap(eid[6], eid[7]);

    F2E[nf0] = { eid[0], ne0   , ne3   , eid[7] };
    F2E[nf1] = { eid[1], eid[2], ne1   , ne0    };
    F2E[nf2] = { ne1   , eid[3], eid[4], ne2    };
    F2E[nf3] = { ne3   , ne2   , eid[5], eid[6] };

    // External E2F
    for (Int i = 0; i < 4; i++) {
        Int fchild = nf + i;
        const auto &pe = F2E[fchild];
        for (Int j = 0; j < 4; j++) {
            Int edge = pe[j];
            if      (E2F[edge][0] == quad) E2F[edge][0] = fchild;
            else if (E2F[edge][1] == quad) E2F[edge][1] = fchild;
        }
    }

    // Internal E2F
    E2F[ne0][0] = nf0; E2F[ne0][1] = nf1;
    E2F[ne1][0] = nf1; E2F[ne1][1] = nf2;
    E2F[ne2][0] = nf2; E2F[ne2][1] = nf3;
    E2F[ne3][0] = nf3; E2F[ne3][1] = nf0;

    // Disable quad and enable its children
    factive.erase(quad);
    factive.insert(nf0);
    factive.insert(nf1);
    factive.insert(nf2);
    factive.insert(nf3);

    // Set quad children pointers
    fchildren[quad] = { nf0, nf1, nf2, nf3 };
    flevel[nf0] = flevel[nf1] = flevel[nf2] = flevel[nf3] = flevel[quad] + 1;

    // Enable internal edges
    eactive.insert(ne0);
    eactive.insert(ne1);
    eactive.insert(ne2);
    eactive.insert(ne3);
    elevel[ne0] = elevel[ne1] = elevel[ne2] = elevel[ne3] = flevel[quad] + 1;

    l2gp[np] = M_np++;
    
    l2ge[ne0] = M_ne++;
    l2ge[ne1] = M_ne++;
    l2ge[ne2] = M_ne++;
    l2ge[ne3] = M_ne++;

    l2gf[nf0] = M_nf++;
    l2gf[nf1] = M_nf++;
    l2gf[nf2] = M_nf++;
    l2gf[nf3] = M_nf++;

    np += 1;
    ne += 4;
    nf += 4;
}

void Smesh::refine_tri(Int tri)
{
    printf(SF_D_ "\n", tri);
    assert(0);
}

void Smesh::get_leaves(Int face, std::vector<Int> &leaves) const
{
    if (face_is_active(face)) {
        leaves.push_back(face);
        return;
    }

    for (Int child : fchildren.at(face)) get_leaves(child, leaves);
}
