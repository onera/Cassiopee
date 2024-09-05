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

#include <cassert>
#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <map>
#include <algorithm>
#include <queue>
#include <stack>

bool Smesh::ccw_oriented(E_Int face)
{
    E_Float sum = 0;

    const auto &pn = F[face];

    for (size_t i = 0; i < pn.size(); i++) {
        E_Int p = pn[i];
        E_Int q = pn[(i+1)%pn.size()];
        sum += (X[q]-X[p])*(Y[q]+Y[p]);
    }

    E_Int sign = Sign(sum);
    assert(sign != 0);

    return sign < 0;
}

Smesh::Smesh(const IMesh &M)
{
    F.resize(M.patch.size());
    //assert(M.patch.size() == M.skin.size());

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

    assert(np == g2lp.size());
    assert(np == l2gp.size());

    // Get the points
    //np = g2lp.size();
    X.resize(np);
    Y.resize(np);
    Z.resize(np);

    for (const auto &pids : g2lp) {
        X[pids.second] = M.X[pids.first];
        Y[pids.second] = M.Y[pids.first];
        Z[pids.second] = M.Z[pids.first];
    }

    make_edges();
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

    assert(ne == E.size());

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


    // Check
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

#define TRI 5
#define QUAD 9

Smesh::Smesh(const char *fname)
{
    FILE *fh = fopen(fname, "r");
    assert(fh);

    char buf[256];

    fgets(buf, 256, fh);
    fgets(buf, 256, fh);
    char *next = strtok(buf, " ");
    next = strtok(NULL, "\n");
   
    char *bad_ptr = NULL;
    nf = strtod(next, &bad_ptr);
    assert(*bad_ptr == '\0'); 
    printf("Faces: " SF_D_ "\n", nf);

    F.resize(nf);

    E_Int ret, type, dummy;

    for (E_Int i = 0; i < nf; i++) {
        auto &cn = F[i];
        ret = fscanf(fh, SF_D_ " ", &type);
        if (ret != 1) abort();
        if (type == TRI) {
            cn.resize(3);
            ret = fscanf(fh, SF_D_ SF_D_ SF_D_ SF_D_ "\n", 
                &cn[0], &cn[1], &cn[2], &dummy);
            if (ret != 4) abort();
        } else if (type == QUAD) {
            cn.resize(4);
            ret = fscanf(fh, SF_D_ SF_D_ SF_D_ SF_D_ SF_D_ "\n",
                &cn[0], &cn[1], &cn[2], &cn[3], &dummy);
            if (ret != 5) abort();
        } else {
            assert(0);
        }
    }

    fgets(buf, 256, fh);
    next = strtok(buf, " ");
    next = strtok(NULL, "\n");
   
    np = strtod(next, &bad_ptr);
    assert(*bad_ptr == '\0'); 
    printf("points: " SF_D_ "\n", np);

    X.resize(np);
    Y.resize(np);

    for (E_Int i = 0; i < np; i++) {
        ret = fscanf(fh, "%lf %lf " SF_D_ "\n", &X[i], &Y[i], &dummy);
        if (ret != 3) abort();
    }

    fclose(fh);

    make_edges();

    init_adaptation_data();
}

// Refine M0/M1 as long as one of its faces contains another face from M1/M0

#define NON_DEGEN 0
#define ON_VERTEX 1
#define ON_EDGE 2

bool Smesh::faces_are_dups(E_Int face, E_Int mface, const Smesh &M)
{
    const auto &pn = F[face];
    const auto &pnm = M.F[mface];

    assert(face_is_quad(face) || face_is_tri(face));
    assert(M.face_is_quad(mface) || M.face_is_tri(mface));

    if (pn.size() != pnm.size()) return false;

    E_Int mfound[4] = { 0, 0, 0, 0 };

    for (size_t i = 0; i < pn.size(); i++) {
        E_Int p = pn[i];
        for (size_t j = 0; j < pnm.size(); j++) {
            E_Int pm = pnm[j];
            if (cmp_points(X[p], Y[p], Z[p], M.X[pm], M.Y[pm], M.Z[pm]) == 0) {
                assert(mfound[i] == 0);
                mfound[i] = 1;
                break;
            }
        }

        if (mfound[i] == 0) return false;
    }

    return true;
}

std::vector<E_Int> Smesh::prepare_for_refinement(const std::vector<E_Int> &ref_data)
{
    std::vector<E_Int> ref_faces;
    
    for (size_t i = 0; i < ref_data.size(); i++) {
        if (ref_data[i] > 0) ref_faces.push_back(i);
    }

    // Refine the lower-level faces first
    std::sort(ref_faces.begin(), ref_faces.end(),
        [&] (E_Int i, E_Int j) { return flevel[i] < flevel[j]; });
    
    // Resize data structures
    resize_point_data(ref_faces.size());
    resize_edge_data(ref_faces.size());
    resize_face_data(ref_faces.size());

    return ref_faces;
}

void Smesh::resize_point_data(size_t nref_faces)
{
    size_t nnew_points = np + nref_faces * 5;
    X.resize(nnew_points);
    Y.resize(nnew_points);
}

void Smesh::resize_edge_data(size_t nref_faces)
{
    size_t nnew_edges = ne + nref_faces * 12;
    E.resize(nnew_edges, {-1, -1});
    E2F.resize(nnew_edges, {-1, -1});
    elevel.resize(nnew_edges, -1);
}

void Smesh::resize_face_data(size_t nref_faces)
{
    size_t nnew_faces = nf + nref_faces * 4;
    F.resize(nnew_faces);
    F2E.resize(nnew_faces);
    flevel.resize(nnew_faces, -1);
}

void Smesh::refine_faces(const std::vector<E_Int> &ref_faces)
{
    for (E_Int ref_face : ref_faces) {

        if (face_is_tri(ref_face)) refine_tri(ref_face);
        else refine_quad(ref_face);
    }
}

std::vector<E_Int> Smesh::smooth_ref_data(std::map<E_Int, std::vector<E_Int>> &sensor)
{
    std::vector<E_Int> ref_data(nf, 0);
    std::stack<E_Int> stk;

    for (const auto &fdata : sensor) {
        assert(face_is_active(fdata.first));
        ref_data[fdata.first] = 1;
        stk.push(fdata.first);
    }

    while (!stk.empty()) {
        E_Int face = stk.top();
        stk.pop();

        E_Int face_incr = ref_data[face] + flevel[face];

        auto neis = get_active_neighbours(face);

        for (auto nei : neis) {
            assert(nei != -1);

            E_Int nei_incr = ref_data[nei] + flevel[nei];

            E_Int diff = abs(face_incr - nei_incr);

            if (diff <= 1) continue;

            E_Int fid = face_incr > nei_incr ? nei : face;

            ref_data[fid] += 1;

            stk.push(fid);
        }
    }

    return ref_data;
}

std::vector<E_Int> Smesh::get_active_neighbours(E_Int face)
{
    std::vector<E_Int> neis;

    for (E_Int edge : F2E[face]) {
        assert(face == E2F[edge][0] || face == E2F[edge][1]);

        // For now, an edge is enabled if it has no children.
        // This is not necessarily true if unrefinement is considered.
        auto it = echildren.find(edge);

        if (it == echildren.end()) {
            assert(edge_is_active(edge));

            E_Int nei = get_neighbour(face, edge);

            if (nei != -1) neis.push_back(nei);
        } else {
            // Edge was refined, get the neighbours of its children
            assert(!edge_is_active(edge));

            const auto &children = it->second;

            for (E_Int child : children)  {
                E_Int nei = get_neighbour(face, child);

                if (nei != -1) neis.push_back(nei);
            }
        }
    }

    return neis;
}

bool Smesh::face_contains_Mface(E_Int face, E_Int mface, const Smesh &M) const
{
    // face containes mface iff it contains all its points
    assert(M.face_is_active(mface));
    const auto &cn = M.F[mface];
    for (E_Int p : cn) {
        if (face_contains_point(face, M.X[p], M.Y[p], M.Z[p]) == -1)
            return false;
    }
    return true;
}

E_Int Smesh::face_contains_point(E_Int face, E_Float x, E_Float y, E_Float z) const
{
    const auto &cn = F[face];

    E_Int hit, a, b, c;

    if (face_is_quad(face)) {

        // First triangle 

        a = cn[0], b = cn[1], c = cn[2];
        
        hit = Triangle::is_point_inside(x, y, z, X[a], Y[a], Z[a],
            X[b], Y[b], Z[b], X[c], Y[c], Z[c]);
        
        if (hit) return 0;

        // Second triangle
        a = cn[0], b = cn[2], c = cn[3];
        
        hit = Triangle::is_point_inside(x, y, z, X[a], Y[a], Z[a],
            X[b], Y[b], Z[b], X[c], Y[c], Z[c]);
        
        if (hit) return 1;
    } else {

        assert(face_is_tri(face));

        a = cn[0], b = cn[1], c = cn[2];
        
        hit = Triangle::is_point_inside(x, y, z, X[a], Y[a], Z[a],
            X[b], Y[b], Z[b], X[c], Y[c], Z[c]);
        
        if (hit) return 0;
    }

    return -1;
}

/*
std::vector<pointFace> Smesh::locate(E_Float x, E_Float y, E_Float z) const
{
    E_Int a, b, c;
    E_Int hit = 0;
    std::vector<pointFace> HITS;

    for (E_Int face : factive) {
        const auto &cn = F[face];
        if (cn.size() == 4) {
            // First triangle
            a = cn[0], b = cn[1], c = cn[2];
            hit = Triangle::is_point_inside(x, y, z, X[a], Y[a], Z[a],
                X[b], Y[b], Z[b], X[c], Y[c], Z[c]);
            if (hit) {
                HITS.push_back(pointFace(face, 0));
            } else {
                // Second triangle
                a = cn[0], b = cn[2], c = cn[3];
                hit = Triangle::is_point_inside(x, y, z, X[a], Y[a], Z[a],
                    X[b], Y[b], Z[b], X[c], Y[c], Z[c]);
                if (hit) {
                    HITS.push_back(pointFace(face, 1));
                }
            }
        } else {
            a = cn[0], b = cn[1], c = cn[2];
            hit = Triangle::is_point_inside(x, y, z, X[a], Y[a], Z[a],
                X[b], Y[b], Z[b], X[c], Y[c], Z[c]);
            if (hit) {
                HITS.push_back(pointFace(face, 0));
            }
        }
    }

    return HITS;
}
*/

void Smesh::init_adaptation_data()
{
    flevel.resize(nf, 0);
    elevel.resize(ne, 0);

    for (E_Int i = 0; i < nf; i++) {
        factive.insert(i);
    }

    for (E_Int i = 0; i < ne; i++) {
        eactive.insert(i);
    }
}

void Smesh::make_point_faces()
{
    P2F.clear();
    P2F.resize(np);

    for (E_Int face : factive) {
        const auto &cn = F[face];
        for (auto p : cn)
            P2F[p].push_back(face);
    }
}

void Smesh::make_point_faces_all()
{
    assert(P2F.empty());
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
}

Smesh Smesh::extract_conformized()
{
    // Keep all the points
    std::vector<E_Float> new_X(np), new_Y(np);
    for (E_Int i = 0; i < np; i++) {
        new_X[i] = X[i];
        new_Y[i] = Y[i];
    }

    // Keep the active edges and faces
    std::map<E_Int, E_Int> new_eids;
    std::map<E_Int, E_Int> new_fids;

    E_Int new_ne = 0, new_nf = 0;

    for (E_Int edge : eactive) new_eids[edge] = new_ne++;
    for (E_Int face : factive) new_fids[face] = new_nf++;

    // E
    std::vector<o_edge> new_E(new_ne, {-1, -1});
    
    for (E_Int edge : eactive) {
        E_Int new_eid = new_eids[edge];

        // TODO: (p,q) -> new_pids(p,q)
        new_E[new_eid].p = E[edge].p;
        new_E[new_eid].q = E[edge].q;
    }

    // E2F
    std::vector<std::array<E_Int, 2>> new_E2F(new_ne, {-1, -1});

    for (E_Int edge : eactive) {
        E_Int new_eid = new_eids[edge];

        E_Int left = E2F[edge][0];
        E_Int right = E2F[edge][1];

        new_E2F[new_eid][0] = new_fids[left];
        if (right != -1) new_E2F[new_eid][1] = new_fids[right];
    }

    // F2E
    std::vector<std::vector<E_Int>> new_F2E(new_nf);

    for (E_Int face : factive) {
        E_Int new_fid = new_fids[face];

        const auto &old_pe = F2E[face];

        auto &new_pe = new_F2E[new_fid];

        for (E_Int e : old_pe) {
            if (edge_is_active(e)) {
                new_pe.push_back(new_eids[e]);
            } else {
                const auto &children = echildren[e];
                assert(children.size() == 2);

                for (E_Int child : children) new_pe.push_back(new_eids[child]);
            }
        }
    }

    // Check consistent E2F and F2E
    /*
    for (E_Int i = 0; i < new_nf; i++) {
        const auto &pe = new_F2E[i];

        for (E_Int e : pe) {
            assert(i == new_E2F[e][0] || i == new_E2F[e][1]);
        }
    }
    */

    // F
    std::vector<std::vector<E_Int>> new_F(new_nf);
    
    for (E_Int face = 0; face < new_nf; face++) {
        const auto &pe = new_F2E[face];

        auto &pn = new_F[face];

        for (E_Int e : pe) {
            assert(face == new_E2F[e][0] || face == new_E2F[e][1]);
            if (new_E2F[e][0] == face) pn.push_back(new_E[e].p);
            else pn.push_back(new_E[e].q);
        }
    }

    // Build and return
    Smesh ret;
    ret.X = new_X;
    ret.Y = new_Y;
    ret.E = new_E;
    ret.E2F = new_E2F;
    ret.F = new_F;
    ret.F2E = new_F2E;
    ret.np = np;
    ret.ne = new_ne;
    ret.nf = new_nf;
    //ret.l2gp = l2gp;
    //ret.factive = factive;

    return ret;
}

void Smesh::write_faces(const char *fname, const std::vector<E_Int> &faces) const
{
    std::map<E_Int, E_Int> pmap;

    E_Int np = 0;

    for (size_t i = 0; i < faces.size(); i++) {
        E_Int f = faces[i];
        const auto &cn = F[f];
        for (auto p : cn) {
            if (pmap.find(p) == pmap.end()) {
                pmap[p] = np++;
            }
        }
    }

    std::vector<E_Int> ipmap(pmap.size());
    for (auto &pdata : pmap) ipmap[pdata.second] = pdata.first;

    FILE *fh = fopen(fname, "w");
    assert(fh);

    fprintf(fh, "POINTS\n");
    fprintf(fh, "%d\n", np);
    for (size_t i = 0; i < ipmap.size(); i++) {
        E_Int op = ipmap[i];
        fprintf(fh, "%f %f %f\n", X[op], Y[op], Z[op]);
    }

    fclose(fh);
}


void Smesh::write_ngon(const char *fname)
{
    FILE *fh = fopen(fname, "w");
    assert(fh);

    fprintf(fh, "POINTS\n");
    fprintf(fh, SF_D_ "\n", np);
    for (E_Int i = 0; i < np; i++) {
        fprintf(fh, "%f %f %f\n", X[i], Y[i], Z[i]);
    }

    fprintf(fh, "INDPG\n");
    fprintf(fh, SF_D_ "\n", ne+1);
    size_t sizeNGon = 0;
    fprintf(fh, "0 ");
    for (E_Int i = 0; i < ne; i++) {
        sizeNGon += 2;
        fprintf(fh, "%zu ", sizeNGon);
    }
    fprintf(fh, "\n");

    fprintf(fh, "NGON\n");
    fprintf(fh, "%zu\n", sizeNGon);
    for (E_Int i = 0; i < ne; i++) {
        fprintf(fh, SF_D_  " " SF_D_  " ", E[i].p, E[i].q);
    }
    fprintf(fh, "\n");

    fprintf(fh, "INDPH\n");
    fprintf(fh, SF_D_ "\n", nf+1);
    size_t sizeNFace = 0;
    fprintf(fh, "0 ");
    for (E_Int i = 0; i < nf; i++) {
        sizeNFace += F2E[i].size();
        fprintf(fh, "%zu ", sizeNFace);
    }
    fprintf(fh, "\n");

    fprintf(fh, "NFACE\n");
    fprintf(fh, "%zu\n", sizeNFace);
    for (E_Int i = 0; i < nf; i++) {
        for (E_Int e : F2E[i]) fprintf(fh, SF_D_ " ", e);
    }
    fprintf(fh, "\n");

    fclose(fh);
}

void Smesh::make_fnormals()
{
    fnormals.clear();
    fnormals.resize(3*nf, 0);
    
    for (E_Int fid = 0; fid < nf; fid++) {
        const auto &pn = F[fid];

        E_Float o[3] = {0, 0, 0};
        for (E_Int p : pn) {
            o[0] += X[p];
            o[1] += Y[p];
            o[2] += Z[p];
        }
        for (int i = 0; i < 3; i++) o[i] /= pn.size();

        E_Int a = pn[0];
        E_Int b = pn[1];

        E_Float oa[3] = {X[a]-o[0], Y[a]-o[1], Z[a]-o[2]};
        E_Float ob[3] = {X[b]-o[0], Y[b]-o[1], Z[b]-o[2]};

        E_Float *N = &fnormals[3*fid];
        K_MATH::cross(oa, ob, N);
        E_Float NORM = K_MATH::norm(N, 3);
        assert(Sign(NORM) != 0);
        for (E_Int i = 0; i < 3; i++) N[i] /= NORM;
    }
}

void Smesh::make_pnormals()
{
    make_fnormals();

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

AABB Smesh::AABB_face(const std::vector<E_Int> &pn) const
{
    AABB ret;
    for (E_Int p : pn) {
        if (X[p] > ret.xmax) ret.xmax = X[p];
        if (X[p] < ret.xmin) ret.xmin = X[p];
        if (Y[p] > ret.ymax) ret.ymax = Y[p];
        if (Y[p] < ret.ymin) ret.ymin = Y[p];
        if (Z[p] > ret.zmax) ret.zmax = Z[p];
        if (Z[p] < ret.zmin) ret.zmin = Z[p];
    }
    return ret;
}

void Smesh::make_bbox()
{
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
}

void Smesh::hash_faces()
{
    NX = 100;
    NY = 100;
    NZ = 10;

    HX = (xmax - xmin) / NX;
    HY = (ymax - ymin) / NY;
    HZ = (zmax - zmin) / NZ;

    NXY = NX * NY;

    fmap.clear();

    for (E_Int fid = 0; fid < nf; fid++) {
        const auto &pn = F[fid];

        AABB bbox = AABB_face(pn);

        E_Int imin = floor((bbox.xmin - xmin) / HX);
        E_Int imax = floor((bbox.xmax - xmin) / HX);
        E_Int jmin = floor((bbox.ymin - ymin) / HY);
        E_Int jmax = floor((bbox.ymax - ymin) / HY);
        E_Int kmin = floor((bbox.zmin - zmin) / HZ);
        E_Int kmax = floor((bbox.zmax - zmin) / HZ);

        for (E_Int k = kmin; k < kmax+1; k++) {
            for (E_Int j = jmin; j < jmax+1; j++) {
                for (E_Int i = imin; i < imax+1; i++) {
                    E_Int voxel = i + NX*j + NXY*k;
                    fmap[voxel].push_back(fid);
                }
            }
        }
    }
}