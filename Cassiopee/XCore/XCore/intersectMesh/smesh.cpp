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

Smesh::Smesh()
{}

Smesh::Smesh(const IMesh &M, bool is_planar)
{
    NEAR_VERTEX_TOL = M.NEAR_VERTEX_TOL;
    NEAR_EDGE_TOL = M.NEAR_EDGE_TOL;

    F.resize(M.patch.size());

    nf = 0;
    np = 0;

    // Get the faces
    for (E_Int gf : M.patch) {
        g2lf[gf] = nf;
        l2gf[nf] = gf;

        auto &face = F[nf];

        for (E_Int gp : M.F[gf]) {
            auto it = g2lp.find(gp);

            if (it == g2lp.end()) {
                g2lp[gp] = np;
                l2gp[np] = gp;
                face.push_back(np);
                np++;
            } else {
                face.push_back(it->second);
            }
        }

        nf++;
    }

    assert((size_t)np == g2lp.size());
    assert((size_t)np == l2gp.size());

    // Get the points
    X.resize(np);
    Y.resize(np);
    Z.resize(np);

    for (const auto &pids : g2lp) {
        X[pids.second] = M.X[pids.first];
        Y[pids.second] = M.Y[pids.first];
        Z[pids.second] = M.Z[pids.first];
    }

    make_edges(is_planar);
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

void Smesh::make_edges(bool is_planar)
{
    // Make the edges
    F2E.resize(F.size());
    std::map<o_edge, E_Int, o_edge_cmp> edges;

    ne = 0;

    for (E_Int i = 0; i < nf; i++) {
        auto &face = F[i];
        for (size_t j = 0; j < face.size(); j++) {
            E_Int p = face[j];
            E_Int q = face[(j+1)%face.size()];
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
    E2F.resize(ne, {-1,-1});

    std::vector<E_Int> count(ne, 0);

    for (E_Int i = 0; i < nf; i++) {
        const auto &face = F2E[i];
        const auto &pn = F[i];
        assert(face.size() == pn.size());

        for (size_t j = 0; j < face.size(); j++) {
            E_Int e = face[j];
            count[e]++;

            if (E2F[e][0] == -1) E2F[e][0] = i;
            else {
                assert(E2F[e][1] == -1);
                E2F[e][1] = i;
            }
        }
    }


    // Check Euler formula for planar graphs
    if (is_planar)
        assert(np - ne + nf + 1 == 2);

    for (E_Int i = 0; i < ne; i++) {
        assert(count[i] == 1 || count[i] == 2);
    }


    for (E_Int i = 0; i < ne; i++) {
        E_Int fi = E2F[i][0];
        E_Int fj = E2F[i][1];

        const auto &face_i = F2E[fi];
        E_Int found = false;
        for (size_t j = 0; j < face_i.size(); j++) {
            if (face_i[j] == (E_Int)i) {
                found = true;
                break;
            }
        }

        if (!found) abort();

        if (fj == -1) continue;

        const auto &face_j = F2E[fj];
        found = false;
        for (size_t j = 0; j < face_j.size(); j++) {
            if (face_j[j] == (E_Int)i) {
                found = true;
                break;
            }
        }

        if (!found) abort();
    }

    // Make faces neighbourhood
    F2F.resize(nf);
    for (size_t i = 0; i < F2E.size(); i++) {
        auto &face = F2E[i];
        auto &neis = F2F[i];
        for (size_t j = 0; j < face.size(); j++) {
            E_Int e = face[j];
            if (E2F[e][0] == (E_Int)i) neis.push_back(E2F[e][1]);
            else if (E2F[e][1] == (E_Int)i) neis.push_back(E2F[e][0]);
            else assert(0);
        }
    }

    // Traverse the face list breadth-first and adjust edges accordingly
    std::vector<E_Int> visited(nf, 0);
    std::queue<E_Int> Q;
    Q.push(0);
    visited[0] = 1;

    while (!Q.empty()) {
        E_Int f = Q.front();
        Q.pop();

        assert(f != -1);

        visited[f] = 1;

        auto &neis = F2F[f];
        auto &edgs = F2E[f];
        auto &face = F[f];

        for (size_t j = 0; j < face.size(); j++) {
            E_Int nei = neis[j];
            
            E_Int p = face[j];
            E_Int q = face[(j+1)%face.size()];
            
            E_Int e = edgs[j];
            
            if (nei == -1) {
                assert(E[e].p == p);
                assert(E[e].q == q);
                continue;
            }

            if (visited[nei]) {
                assert(E2F[e][0] == nei);
                assert(E2F[e][1] == f);
                assert(E[e].p == q);
                assert(E[e].q == p);
                continue;
            }
 
            if (E[e].p != p) {
                assert(visited[nei] == 0);
                assert(E[e].q == p);
                assert(E[e].p == q);
                std::swap(E[e].p, E[e].q);
                E2F[e][0] = f;
                E2F[e][1] = nei;
                Q.push(nei);
            } 
        }
    }

    // Check
    for (E_Int i = 0; i < nf; i++) {
        auto &face = F[i];
        for (size_t j = 0; j < face.size(); j++) {
            E_Int e = F2E[i][j];
            E_Int p = face[j];
            E_Int q = face[(j+1)%face.size()];

            if (E[e].p == p) {
                assert(E[e].q == q);
                assert(E2F[e][0] == (E_Int)i);
            } else if (E[e].q == p) {
                assert(E[e].p == q);
                assert(E2F[e][1] == (E_Int)i);
            } else {
                assert(0);
            }
        }
    }
}

void Smesh::make_fcenters()
{
    fcenters.clear();
    fcenters.resize(3*nf, 0);
    for (E_Int fid = 0; fid < nf; fid++) {
        E_Float *fc = &fcenters[3*fid];
        const auto &pn = F[fid];
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

    for (E_Int face = 0; face < nf; face++) {
        const auto &cn = F[face];
        for (auto p : cn)
            P2F[p].push_back(face);
    }
}

void Smesh::make_point_edges()
{
    P2E.clear();
    P2E.resize(np);

    for (E_Int eid = 0; eid < ne; eid++) {
        P2E[E[eid].p].push_back(eid);
        P2E[E[eid].q].push_back(eid);
    }

    for (E_Int pid = 0; pid < np; pid++) {
        const auto &pe = P2E[pid];
        for (auto eid : pe) {
            const auto &e = E[eid];
            assert(e.p == pid || e.q == pid);
        }
    }

}

void Smesh::make_fnormals()
{
    fnormals.clear();
    fnormals.resize(3*nf, 0);
    
    for (E_Int fid = 0; fid < nf; fid++) {
        const auto &pn = F[fid];
        E_Int a = pn[0], b = pn[1], c = pn[2];
        E_Float v0[3] = {X[b]-X[a], Y[b]-Y[a], Z[b]-Z[a]};
        E_Float v1[3] = {X[c]-X[a], Y[c]-Y[a], Z[c]-Z[a]};
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

    for (E_Int pid = 0; pid < np; pid++) {
        const auto &faces = P2F[pid];

        E_Float *N = &pnormals[3*pid];

        for (E_Int fid : faces) {
            E_Float *fN = &fnormals[3*fid];
            for (E_Int i = 0; i < 3; i++) N[i] += fN[i];
        }

        E_Float NORM = K_MATH::norm(N, 3);

        for (E_Int i = 0; i < 3; i++) N[i] /= NORM;
    }
}

void Smesh::make_bbox()
{
    NX = 100;
    NY = 100;
    NZ = 100;
    NXY = NX * NY;
    NXYZ = NXY * NZ;

    xmin = ymin = zmin = std::numeric_limits<E_Float>::max();
    xmax = ymax = zmax = std::numeric_limits<E_Float>::min();

    for (E_Int i = 0; i < np; i++) {
        if (X[i] < xmin) xmin = X[i];
        if (Y[i] < ymin) ymin = Y[i];
        if (Z[i] < zmin) zmin = Z[i];
        if (X[i] > xmax) xmax = X[i];
        if (Y[i] > ymax) ymax = Y[i];
        if (Z[i] > zmax) zmax = Z[i];
    }

    E_Float dx = xmax - xmin;
    E_Float dy = ymax - ymin;
    E_Float dz = zmax - zmin;

    xmin = xmin - dx*0.01;
    ymin = ymin - dy*0.01;
    zmin = zmin - dz*0.01;
    xmax = xmax + dx*0.01;
    ymax = ymax + dy*0.01;
    zmax = zmax + dz*0.01;

    HX = (xmax - xmin) / NX;
    HY = (ymax - ymin) / NY;
    HZ = (zmax - zmin) / NZ;
}

void Smesh::hash_faces()
{
    bin_faces.clear();
    bin_faces.resize(NXYZ);

    for (E_Int fid = 0; fid < nf; fid++) {
        const auto &pn = F[fid];

        E_Int Imin, Jmin, Kmin;
        E_Int Imax, Jmax, Kmax;

        Imin = Jmin = Kmin = NXYZ;
        Imax = Jmax = Kmax = -1;

        for (E_Int p : pn) {
            E_Float x = X[p];
            E_Float y = Y[p];
            E_Float z = Z[p];

            E_Int I = floor((x - xmin) / HX);
            E_Int J = floor((y - ymin) / HY);
            E_Int K = floor((z - zmin) / HZ);

            if (I < Imin) Imin = I;
            if (J < Jmin) Jmin = J;
            if (K < Kmin) Kmin = K;
            if (I > Imax) Imax = I;
            if (J > Jmax) Jmax = J;
            if (K > Kmax) Kmax = K;
        }

        for (E_Int I = Imin; I <= Imax; I++) {
            for (E_Int J = Jmin; J <= Jmax; J++) {
                for (E_Int K = Kmin; K <= Kmax; K++) {
                    E_Int voxel = get_voxel(I, J, K);
                    assert(voxel >= 0);
                    assert(voxel < NXYZ);
                    bin_faces[voxel].push_back(fid);
                }
            }
        }
    }
}

void Smesh::compute_min_distance_between_points()
{
    min_pdist_squared = EFLOATMAX;
    E_Int ndists = 0;

    for (E_Int i = 0; i < np; i++) {
        E_Float xi = X[i];
        E_Float yi = Y[i];
        E_Float zi = Z[i];
        for (E_Int j = i+1; j < np; j++) {
            E_Float xj = X[j];
            E_Float yj = Y[j];
            E_Float zj = Z[j];

            E_Float dx = xj-xi;
            E_Float dy = yj-yi;
            E_Float dz = zj-zi;

            E_Float dist = dx*dx + dy*dy + dz*dz;

            if (dist < min_pdist_squared) min_pdist_squared = dist;

            ndists++;
        }
    }

    assert(ndists == np*(np-1)/2);
}

