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
#include <cstdio>
#include <cassert>
#include <cstring>
#include <cstdlib>
#include <limits>
#include <map>
#include <set>
#include <cmath>
#include <algorithm>
#include <ctime>

#include "mesh.h"
#include "triangle.h"
#include "primitives.h"
#include "ray.h"
#include "io.h"
#include "common/Karray.h"

void IMesh::triangulate(const std::vector<E_Int> &fids)
{
    owner.clear();
    neigh.clear();
    owner.resize(nf + fids.size(), -1);
    neigh.resize(nf + fids.size(), -1);
    
    for (E_Int cid = 0; cid < nc; cid++) {
        const auto &pf = C[cid];
        for (E_Int fid : pf) {
            if (owner[fid] == -1) owner[fid] = cid;
            else neigh[fid] = cid;
        }
    }

    for (auto fid : fids) {
        const auto &pn = F[fid];

        assert(pn.size() == 4);

        E_Int nodes[4] = {pn[0], pn[1], pn[2], pn[3]};

        F.push_back({nodes[2], nodes[3], nodes[0]});
        F[fid] = {nodes[0], nodes[1], nodes[2]};

        E_Int own = owner[fid];
        auto &pf = C[own];
        pf.push_back(nf);

        E_Int nei = neigh[fid];
        if (nei != -1) {
            auto &pf = C[nei];
            pf.push_back(nf);
        }

        owner[nf] = own;
        neigh[nf] = nei;

        nf++;
    }
}

void IMesh::triangulate_face_set(bool propagate)
{
    make_skin();
    
    if (propagate) {
        E_Float xmin, ymin, zmin, xmax, ymax, zmax;
        xmin = ymin = zmin = EFLOATMAX;
        xmax = ymax = zmax = EFLOATMIN;

        for (auto fid : faces_to_tri) {
            const auto &pn = F[fid];
            for (auto p : pn) {
                if (X[p] < xmin) xmin = X[p];
                if (Y[p] < ymin) ymin = Y[p];
                if (Z[p] < zmin) zmin = Z[p];
                if (X[p] > xmax) xmax = X[p];
                if (Y[p] > ymax) ymax = Y[p];
                if (Z[p] > zmax) zmax = Z[p];
            }
        }


        faces_to_tri.clear();

        for (auto fid : skin) {
            const auto &pn = F[fid];
            for (auto p : pn) {
                if (X[p] >= xmin && X[p] <= xmax &&
                    Y[p] >= ymin && Y[p] <= ymax &&
                    Z[p] >= zmin && Z[p] <= zmax) {
                    faces_to_tri.insert(fid);
                    break;
                }
            }
        }
    }

    E_Int NF = nf;

    E_Int face_incr = (E_Int)faces_to_tri.size();

    F.resize(NF + face_incr);

    patch.clear();

    for (E_Int fid : faces_to_tri) {

        auto &pn = F[fid];

        patch.insert(fid);

        if (pn.size() == 3) continue;

        assert(pn.size() == 4);

        std::vector<E_Int> tri0(3), tri1(3);

        tri0[0] = pn[0], tri0[1] = pn[1], tri0[2] = pn[2];

        tri1[0] = pn[2], tri1[1] = pn[3], tri1[2] = pn[0];

        pn = tri0;

        F[NF] = tri1;

        E_Int own = owner[fid];

        assert(own != -1);

        auto &pf = C[own];

        pf.push_back(NF);

        E_Int nei = neigh[fid];

        assert(nei == -1);

        skin.push_back(NF);

        // Update patch

        patch.insert(NF);

        NF++;
    }

    nf = NF;
    F.resize(NF);
}

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

void IMesh::make_edges()
{
    F2E.clear();
    F2E.resize(F.size());
    std::map<o_edge, E_Int, o_edge_cmp> edges;

    E.clear();
    assert(E.empty());
    ne = 0;

    for (E_Int fid = 0; fid < nf; fid++) {
        auto &pn = F[fid];
        for (size_t j = 0; j < pn.size(); j++) {
            E_Int p = pn[j];
            E_Int q = pn[(j+1)%pn.size()];
            o_edge EDGE(p, q);
            auto it = edges.find(EDGE);
            if (it == edges.end()) {
                F2E[fid].push_back(ne);
                edges[EDGE] = ne;
                E.push_back(EDGE);
                ne++;
            } else {
                F2E[fid].push_back(it->second);
            }
        }
    }

    assert((size_t)ne == E.size());

    E2F.clear();
    E2F.resize(ne);

    for (E_Int fid = 0; fid < nf; fid++) {
        const auto &pe = F2E[fid];

        for (E_Int eid : pe) {
            E2F[eid].push_back(fid);
        }
    }
}

E_Float IMesh::get_min_edge_length() const
{
    E_Float ret = EFLOATMAX;
    for (const auto &e : E) {
        E_Int p = e.p;
        E_Int q = e.q;
        E_Float d = (X[p]-X[q])*(X[p]-X[q]) +
                    (Y[p]-Y[q])*(Y[p]-Y[q]) +
                    (Z[p]-Z[q])*(Z[p]-Z[q]);
        d = sqrt(d);
        if (d < ret) ret = d;
    }
    return ret;
}

IMesh::IMesh()
{}

IMesh::IMesh(const Karray &karray)
{
    np = karray.npoints();
    nf = karray.nfaces();
    nc = karray.ncells();

    X.resize(np);
    Y.resize(np);
    Z.resize(np);
    for (E_Int i = 0; i < np; i++) {
        X[i] = karray.x[i];
        Y[i] = karray.y[i];
        Z[i] = karray.z[i];
    }

    F.reserve(nf);
    for (E_Int fid = 0; fid < nf; fid++) {
        E_Int np = -1;
        E_Int *pn = karray.get_face(fid, np);
        std::vector<E_Int> points(np);
        for (E_Int j = 0; j < np; j++)
            points[j] = pn[j] - 1;
        F.push_back(points);
    }

    C.reserve(nc);
    for (E_Int cid = 0; cid < nc; cid++) {
        E_Int nf = -1;
        E_Int *pf = karray.get_cell(cid, nf);
        std::vector<E_Int> faces(nf);
        for (E_Int j = 0; j < nf; j++)
            faces[j] = pf[j] - 1;
        C.push_back(faces);
    }
}

void IMesh::make_point_faces()
{
    P2F.clear();
    P2F.resize(np);
    
    for (E_Int face : factive) {
        const auto &pn = F[face];
        for (E_Int p : pn) P2F[p].push_back(face);
    }
}

void IMesh::make_bbox()
{
    // Grid
    NX = 100;
    NY = 100;
    NZ = 100;
    NXY = NX * NY;
    NXYZ = NXY * NZ;

    xmin = ymin = zmin = EFLOATMAX;
    xmax = ymax = zmax = EFLOATMIN;

    for (E_Int i = 0; i < np; i++) {
        if (X[i] < xmin) xmin = X[i];
        if (Y[i] < ymin) ymin = Y[i];
        if (Z[i] < zmin) zmin = Z[i];
        if (X[i] > xmax) xmax = X[i];
        if (Y[i] > ymax) ymax = Y[i];
        if (Z[i] > zmax) zmax = Z[i];
    }

    xmin = xmin - (xmax - xmin) * 0.01;
    ymin = ymin - (ymax - ymin) * 0.01;
    zmin = zmin - (zmax - zmin) * 0.01;
    xmax = xmax + (xmax - xmin) * 0.01;
    ymax = ymax + (ymax - ymin) * 0.01;
    zmax = zmax + (zmax - zmin) * 0.01;

    HX = (xmax - xmin) / NX;
    HY = (ymax - ymin) / NY;
    HZ = (zmax - zmin) / NZ;
}

void IMesh::make_skin()
{
    skin.clear();

    owner.clear();
    neigh.clear();

    owner.resize(nf, -1);
    neigh.resize(nf, -1);

    for (E_Int cid = 0; cid < nc; cid++) {
        const auto &pf = C[cid];
        for (auto fid : pf) {
           if (owner[fid] == -1) owner[fid] = cid;
           else neigh[fid] = cid;
        }
    }

    for (E_Int fid = 0; fid < nf; fid++) {
        if (neigh[fid] == -1) skin.push_back(fid);
    }
}

void IMesh::write_face(const char *fname, E_Int fid) const
{
    FILE *fh = fopen(fname, "w");
    assert(fh);

    const auto &pn = F[fid];

    fprintf(fh, "POINTS\n");
    fprintf(fh, "%zu\n", pn.size());

    for (E_Int p : pn) {
        fprintf(fh, "%f %f %f\n", X[p], Y[p], Z[p]);
    }

    fclose(fh);
}

bool IMesh::face_contains_sface(E_Int face, E_Int sface, const IMesh &S) const
{
    // face containes mface iff it contains all its points
    assert(S.face_is_active(sface));
    const auto &cn = S.F[sface];
    for (E_Int p : cn) {
        if (face_contains_point(face, S.X[p], S.Y[p], S.Z[p]) == 0) 
            return false;
    }
    return true;
}

E_Int IMesh::face_contains_point(E_Int face, E_Float x, E_Float y, E_Float z) const
{
    const auto &cn = F[face];

    E_Float o[3] = {};

    for (E_Int p : cn) {
        o[0] += X[p];
        o[1] += Y[p];
        o[2] += Z[p];
    }
    for (E_Int i = 0; i < 3; i++) o[i] /= cn.size();

    E_Int hit = 0;

    for (size_t i = 0; i < cn.size(); i++) {
        E_Int a = cn[i];
        E_Int b = cn[(i+1)%cn.size()];

        hit = Triangle::is_point_inside(x, y, z,
            X[a], Y[a], Z[a],
            X[b], Y[b], Z[b],
            o[0], o[1], o[2]);

        if (hit) break;
    }

    return hit;
}

void IMesh::extract_edge_points(E_Int a, E_Int b, std::list<E_Int> &points)
{
    E_Int ref = 0;

    points.clear();
    points.push_back(a);
    points.push_back(b);

    do {
        ref = 0;

        assert(*std::prev(points.end()) == b);

        for (auto it = points.begin(); it != std::prev(points.end()); it++) {
            E_Int a = *it;
            E_Int b = *std::next(it);

            UEdge e(a, b);

            auto search = ecenter.find(e);

            if (search != ecenter.end()) {
                points.insert(std::next(it), search->second);
                ref = 1;
            }
        }
    } while (ref);
}

void IMesh::get_fleaves(E_Int face, std::vector<E_Int> &fleaves)
{
    if (face_is_active(face)) {
        fleaves.push_back(face);
        return;
    }

    for (E_Int child : fchildren.at(face)) get_fleaves(child, fleaves);
}
