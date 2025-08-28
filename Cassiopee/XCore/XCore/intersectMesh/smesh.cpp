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
#include "primitives.h"
#include "triangle.h"
#include "mesh.h"
#include "io.h"
#include "point.h"

#include <cassert>
#include <cstring>
#include <cstdlib>
#include <map>
#include <algorithm>
#include <queue>
#include <stack>

void Smesh::get_edge_centers(E_Int p, E_Int q, std::vector<E_Int> &edge_centers)
{
    u_edge e(p, q);
    auto it = ecenter.find(e);
    if (it == ecenter.end()) 
    {
        edge_centers.push_back(p);
        return;
    }

    E_Int ec = it->second;
    get_edge_centers(p, ec, edge_centers);
    get_edge_centers(ec, q, edge_centers);
}

void Smesh::clear_conformal_data()
{
    //Fc.clear();
    F2E.clear();
    F2F.clear();

    P2F.clear();
    P2E.clear();

    ne = 0;
    E.clear();
    E2F.clear();
}

void Smesh::conformize()
{
    clear_conformal_data();

    std::vector<std::vector<E_Int>> new_Fc(nf);
    
    for (E_Int fid = 0; fid < nf; fid++) 
    {
        const auto &pn = Fc[fid];

        std::vector<E_Int> new_pn;

        for (size_t i = 0; i < pn.size(); i++) 
        {
            E_Int p = pn[i];
            E_Int q = pn[(i+1)%pn.size()];

            std::vector<E_Int> local;
            get_edge_centers(p, q, local);
            for (E_Int p : local) new_pn.push_back(p);
            assert(local[0] == p);
        }

        new_Fc[fid] = new_pn;
    }

    Fc = new_Fc;

    make_edges();
    make_point_faces();
    make_point_edges();
}

Smesh::Smesh()
{}

Smesh Smesh::Smesh_from_mesh_skin(const IMesh &M,
    const std::vector<E_Int> &skin, bool check_Euler)
{
    return Smesh(M, skin, check_Euler);
}

Smesh Smesh::Smesh_from_point_tags(const IMesh &M, const E_Float *ptag,
    bool check_Euler)
{
    std::vector<E_Int> fids;
    for (E_Int fid = 0; fid < M.nf; fid++) 
    {
        const auto &pn = M.F[fid];
        bool keep = true;
        for (auto p : pn) 
        {
            if (ptag[p] != 1.0) 
            {
                keep = false;
                break;
            }
        }
        if (keep) 
        {
            assert(pn.size() == 3);
            fids.push_back(fid);
        }
    }

    return Smesh(M, fids, check_Euler);
}

Smesh Smesh::Smesh_from_mesh_patch(const IMesh &M, bool check_Euler)
{
    std::vector<E_Int> fids;

    for (E_Int fid : M.patch) {
        fids.push_back(fid);
    }

    return Smesh(M, fids, check_Euler);
}

Smesh Smesh::Smesh_from_tagged_faces(const IMesh &M, bool check_Euler)
{
    return Smesh(M, M.ftag, check_Euler);
}

Smesh::Smesh(const IMesh &M, const std::vector<E_Int> &fids, bool check_Euler_)
{
    NEAR_VERTEX_TOL = M.NEAR_VERTEX_TOL;
    NEAR_EDGE_TOL = M.NEAR_EDGE_TOL;

    check_Euler = check_Euler_;

    F.resize(fids.size());

    nf = 0;
    np = 0;

    // Get the faces
    for (E_Int gf : fids) {
        g2lf[gf] = nf;
        l2gf[nf] = gf;

        auto &pn = F[nf];

        for (E_Int gp : M.F[gf]) {
            auto it = g2lp.find(gp);

            if (it == g2lp.end()) {
                g2lp[gp] = np;
                l2gp[np] = gp;
                pn.push_back(np);
                np++;
            } else {
                pn.push_back(it->second);
            }
        }

        nf++;
    }

    assert((size_t)np == g2lp.size());
    assert((size_t)np == l2gp.size());

    // Cache the original number of points and faces
    np_before_adapt = np;
    nf_before_adapt = nf;

    // Get the points
    X.resize(np);
    Y.resize(np);
    Z.resize(np);

    for (const auto &pids : g2lp) {
        X[pids.second] = M.X[pids.first];
        Y[pids.second] = M.Y[pids.first];
        Z[pids.second] = M.Z[pids.first];
    }

    Fc = F;
    make_edges();
    make_point_faces();
    make_point_edges();
    make_bbox();
}

Smesh::Smesh(const Smesh &Mf, const std::set<E_Int> &fids, bool check_Euler_)
{
    NEAR_VERTEX_TOL = Mf.NEAR_VERTEX_TOL;
    NEAR_EDGE_TOL = Mf.NEAR_EDGE_TOL;

    check_Euler = check_Euler_;

    F.resize(fids.size());

    nf = 0;
    np = 0;

    // Get the faces
    for (E_Int gf : fids) {
        g2lf[gf] = nf;
        l2gf[nf] = gf;

        auto &pn = F[nf];

        for (E_Int gp : Mf.Fc[gf]) {
            auto it = g2lp.find(gp);

            if (it == g2lp.end()) {
                g2lp[gp] = np;
                l2gp[np] = gp;
                pn.push_back(np);
                np++;
            } else {
                pn.push_back(it->second);
            }
        }

        nf++;
    }

    assert((size_t)np == g2lp.size());
    assert((size_t)np == l2gp.size());

    // Cache the original number of points and faces
    np_before_adapt = np;
    nf_before_adapt = nf;

    // Get the points
    X.resize(np);
    Y.resize(np);
    Z.resize(np);

    for (const auto &pids : g2lp) {
        X[pids.second] = Mf.X[pids.first];
        Y[pids.second] = Mf.Y[pids.first];
        Z[pids.second] = Mf.Z[pids.first];
    }

    Fc = F;
    make_edges();
    make_point_faces();
    make_point_edges();
    make_bbox();
}

o_edge::o_edge(E_Int P, E_Int Q)
: p(P), q(Q)
{}

struct o_edge_cmp {
    bool operator()(const o_edge &e, const o_edge &f) const
    {
        E_Int e_p = std::min(e.p, e.q);
        E_Int e_q = std::max(e.p, e.q);
        E_Int f_p = std::min(f.p, f.q);
        E_Int f_q = std::max(f.p, f.q);
        return (e_p < f_p) ||
               (e_p == f_p && e_q < f_q);
    }
};

void Smesh::make_edges()
{
    // Make the edges
    F2E.resize(Fc.size());
    std::map<o_edge, E_Int, o_edge_cmp> edges;

    assert(E.empty());
    ne = 0;

    for (E_Int i = 0; i < nf; i++) {
        auto &pn = Fc[i];
        for (size_t j = 0; j < pn.size(); j++) {
            E_Int p = pn[j];
            E_Int q = pn[(j+1)%pn.size()];
            o_edge EDGE(p, q);
            auto it = edges.find(EDGE);
            if (it == edges.end()) {
                F2E[i].push_back(ne);
                edges[EDGE] = ne;
                E.push_back(EDGE);
                ne++;
            } else {
                F2E[i].push_back(it->second);
            }
        }
    }

    assert((size_t)ne == E.size());

    // Make edge_to_face
    assert(E2F.empty());
    E2F.resize(ne, {-1,-1});

    std::vector<E_Int> count(ne, 0);

    for (E_Int i = 0; i < nf; i++) 
    {
        const auto &pe = F2E[i];
        assert(pe.size() == Fc[i].size());

        for (size_t j = 0; j < pe.size(); j++) {
            E_Int e = pe[j];
            count[e]++;

            if (E2F[e][0] == -1) E2F[e][0] = i;
            else 
            {
                assert(E2F[e][1] == -1);
                E2F[e][1] = i;
            }
        }
    }


    // Check Euler formula for planar graphs
    if (check_Euler)
        assert(np - ne + nf + 1 == 2);

    // Check that each edge's count is 1 or 2
    for (E_Int i = 0; i < ne; i++) {
        assert(count[i] == 1 || count[i] == 2);
    }

    // Check consistency between E2F and F2E
    for (E_Int i = 0; i < ne; i++) {
        E_Int fi = E2F[i][0];
        E_Int fj = E2F[i][1];

        const auto &pe_i = F2E[fi];
        E_Int found = false;
        for (size_t j = 0; j < pe_i.size(); j++) {
            if (pe_i[j] ==  i) {
                found = true;
                break;
            }
        }

        if (!found) abort();

        if (fj == -1) continue;

        const auto &pe_j = F2E[fj];
        found = false;
        for (size_t j = 0; j < pe_j.size(); j++) {
            if (pe_j[j] == i) {
                found = true;
                break;
            }
        }

        if (!found) abort();
    }

    // Make faces neighbourhood
    assert(F2F.empty());
    F2F.resize(nf);
    for (E_Int i = 0; i < nf; i++) {
        const auto &pe = F2E[i];
        auto &neis = F2F[i];
        for (size_t j = 0; j < pe.size(); j++) {
            E_Int e = pe[j];
            assert(E2F[e][0] == i || E2F[e][1] == i);
            if (E2F[e][0] == i) neis.push_back(E2F[e][1]);
            else neis.push_back(E2F[e][0]);
        }
    }
}

void Smesh::make_fcenters()
{
    fcenters.clear();
    fcenters.resize(3*nf, 0);
    for (E_Int fid = 0; fid < nf; fid++) {
        E_Float *fc = &fcenters[3*fid];
        const auto &pn = Fc[fid];
        for (E_Int p : pn) {
            fc[0] += X[p];
            fc[1] += Y[p];
            fc[2] += Z[p];
        }
        for (E_Int i = 0; i < 3; i++) fc[i] /= pn.size();
    }
}

void Smesh::make_point_faces()
{
    P2F.clear();
    P2F.resize(np);

    for (E_Int fid = 0; fid < nf; fid++) 
    {
        const auto &pn = Fc[fid];
        for (auto p : pn) P2F[p].push_back(fid);
    }
}

void Smesh::make_point_edges()
{
    P2E.clear();
    P2E.resize(np);

    for (E_Int eid = 0; eid < ne; eid++) 
    {
        P2E[E[eid].p].push_back(eid);
        P2E[E[eid].q].push_back(eid);
    }

#ifndef NDEBUG
    for (E_Int pid = 0; pid < np; pid++) 
    {
        const auto &pe = P2E[pid];
        for (auto eid : pe)
        {
            assert(E[eid].p == pid || E[eid].q == pid);
        }
    }
#endif

}

void Smesh::make_fnormals()
{
    fnormals.clear();
    fnormals.resize(3*nf, 0);
    
    for (E_Int fid = 0; fid < nf; fid++) {
        const auto &pn = Fc[fid];
        const E_Float *fc = &fcenters[3*fid];
        E_Int a = pn[0], b = pn[1];
        E_Float v0[3] = {X[a]-fc[0], Y[a]-fc[1], Z[a]-fc[2]};
        E_Float v1[3] = {X[b]-fc[0], Y[b]-fc[1], Z[b]-fc[2]};
        E_Float *N = &fnormals[3*fid];
        K_MATH::cross(v0, v1, N);
        E_Float NORM = K_MATH::norm(N, 3);
        for (E_Int i = 0; i < 3; i++) N[i] /= NORM;
    }
}

void Smesh::make_pnormals()
{
    pnormals.clear();
    pnormals.resize(3*np, 0);
    
    // Point normals: aggregate of shared faces normals

    assert(fnormals.size() == (size_t)(3*nf));

    for (E_Int pid = 0; pid < np; pid++) {
        const auto &pf = P2F[pid];

        E_Float *N = &pnormals[3*pid];

        for (E_Int fid : pf) {
            E_Float *fN = &fnormals[3*fid];
            for (E_Int i = 0; i < 3; i++) N[i] += fN[i];
        }

        E_Float NORM = K_MATH::norm(N, 3);

        for (E_Int i = 0; i < 3; i++) N[i] /= NORM;
    }
}

void Smesh::clear()
{
    np = ne = nf = 0;
    X.clear();
    Y.clear();
    Z.clear();
    P2F.clear();
    P2E.clear();
    E.clear();
    E2F.clear();
    F.clear();
    F2E.clear();
    F2F.clear();
    g2lp.clear();
    l2gp.clear();
    g2lf.clear();
    l2gf.clear();
}

void Smesh::tag_faces(IMesh &M) const
{
    M.ftag.clear();
    M.ftag.reserve(nf);

    for (E_Int fid = 0; fid < nf; fid++) {
        E_Int gfid = l2gf.at(fid);
        assert(gfid < M.nf);
        M.ftag.push_back(gfid);
    }
}
