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

struct DEdge {
    E_Int p, q;

    DEdge(E_Int P, E_Int Q) : p(P), q(Q) {}

    bool operator<(const DEdge &f) const
    {
        E_Int ep = std::min(p, q);
        E_Int eq = std::max(p, q);
        E_Int fp = std::min(f.p, f.q);
        E_Int fq = std::max(f.p, f.q);
        return (ep < fp) || (ep == fp && eq < fq);
    }
};

void Mesh::make_edges()
{
    std::map<DEdge, E_Int> edges;

    ne = 0;

    F2E.resize(nf);

    for (E_Int i = 0; i < nf; i++) {
        const auto &pn = F[i];
        
        F2E[i].resize(pn.size());

        for (size_t j = 0; j < pn.size(); j++) {
            E_Int p = pn[j];
            E_Int q = pn[(j+1)%pn.size()];
            DEdge e(p, q); 
            auto it = edges.find(e);
            if (it == edges.end()) {
                F2E[i][j] = ne;
                edges[e] = ne;
                ne++;
            } else {
                F2E[i][j] = it->second;
            }
        }
    }

    E.resize(ne);

    for (const auto &edata : edges) {
        E[edata.second][0] = edata.first.p;
        E[edata.second][1] = edata.first.q;
    }
}


std::vector<pointFace> Mesh::locate(E_Float px, E_Float py,
    const std::set<E_Int> &patch) const
{
    E_Int a, b, c;
    E_Int hit;
    std::vector<pointFace> hits;

    for (E_Int face : patch) {
        assert(face_is_active(face));

        const auto &cn = F[face];

        // First triangle

        a = cn[0], b = cn[1], c = cn[2];

        hit = Triangle::ispointInside(px, py, X[a], Y[a], X[b], Y[b], X[c],
            Y[c]);

        if (hit) {
            hits.push_back(pointFace(face, 0));
            continue;
        }
        
        if (face_is_tri(face)) continue;

        // Second triangle
        a = cn[0], b = cn[2], c = cn[3];

        hit = Triangle::ispointInside(px, py, X[a], Y[a], X[b], Y[b], X[c],
            Y[c]);

        if (hit) {
            hits.push_back(pointFace(face, 1));
        }
    }

    return hits;
}

void Mesh::init_adaptation_data()
{
    flevel.resize(nf, 0);

    for (E_Int i = 0; i < nf; i++) factive.insert(i);
}

bool Mesh::faces_are_dups(E_Int mface, E_Int sface, const Mesh &S)
{
    const auto &pnm = F[mface];
    const auto &pns = S.F[sface];

    assert(face_is_quad(mface) || face_is_tri(mface));
    assert(S.face_is_quad(sface) || S.face_is_tri(sface));

    if (pnm.size() != pns.size()) return false;

    E_Int mfound[4] = { 0, 0, 0, 0 };

    for (size_t i = 0; i < pnm.size(); i++) {
        E_Int pm = pnm[i];
        for (size_t j = 0; j < pns.size(); j++) {
            E_Int ps = pns[j];
            if (cmp_points(X[pm], Y[pm], S.X[ps], S.Y[ps]) == 0) {
                assert(mfound[i] == 0);
                mfound[i] = 1;
                break;
            }
        }

        if (mfound[i] == 0) return false;
    }

    return true;
}

Mesh::Mesh()
{}

Mesh::Mesh(K_FLD::FldArrayI &cn, E_Float *x, E_Float *y, E_Float *z, E_Int npts)
{
    np = npts;
    ne = 0;
    nf = cn.getNFaces();
    nc = cn.getNElts();

    X.resize(np);
    Y.resize(np);
    Z.resize(np);
    for (E_Int i = 0; i < np; i++) {
        X[i] = x[i];
        Y[i] = y[i];
        Z[i] = z[i];
    }

    F.reserve(nf);
    E_Int *indPG = cn.getIndPG();
    E_Int *ngon = cn.getNGon();
    for (E_Int i = 0; i < nf; i++) {
        E_Int np = -1;
        E_Int *pn = cn.getFace(i, np, ngon, indPG);
        std::vector<E_Int> points(np);
        for (E_Int j = 0; j < np; j++)
            points[j] = pn[j] - 1;
        F.push_back(points);
    }

    C.reserve(nc);
    E_Int *indPH = cn.getIndPH();
    E_Int *nface = cn.getNFace();
    for (E_Int i = 0; i < nc; i++) {
        E_Int nf = -1;
        E_Int *pf = cn.getElt(i, nf, nface, indPH);
        std::vector<E_Int> faces(nf);
        for (E_Int j = 0; j < nf; j++)
            faces[j] = pf[j] - 1;
        C.push_back(faces);
    }

    make_skin();

    make_bbox();

    hash_skin();

    make_point_faces();

    //make_edges();

    init_adaptation_data();

    srand(time(NULL));
}

void Mesh::make_point_faces()
{
    P2F.clear();
    P2F.resize(np);
    
    for (E_Int face : factive) {
        const auto &pn = F[face];
        for (E_Int p : pn) P2F[p].push_back(face);
    }
}

Mesh::Mesh(const char *fname)
{
    FILE *fh = fopen(fname, "r");
    assert(fh);

    np = ne = nf = nc = 0;

    char buf[256];

    // POINTS
    fgets(buf, 256, fh);
    fgets(buf, 256, fh);
    char *next = strtok(buf, "\n");

    char *bad_ptr = NULL;
    np = strtod(next, &bad_ptr);
    assert(*bad_ptr == '\0');
    printf("points: %d\n", np);

    X.resize(np);
    Y.resize(np);
    Z.resize(np);

    E_Int ret;

    for (E_Int i = 0; i < np; i++) {
        ret = fscanf(fh, "%lf %lf %lf\n", &X[i], &Y[i], &Z[i]);
        assert(ret == 3);
    }

    // FACES
    fgets(buf, 256, fh);
    fgets(buf, 256, fh);
    next = strtok(buf, "\n");

    bad_ptr = NULL;
    nf = strtod(next, &bad_ptr);
    assert(*bad_ptr == '\0');
    printf("Faces: %d\n", nf);

    F.resize(nf);

    for (E_Int i = 0; i < nf; i++) {
        E_Int stride;
        ret = fscanf(fh, "%d ", &stride);
        assert(ret == 1);
        auto &cn = F[i];
        cn.resize(stride);
        for (E_Int j = 0; j < stride-1; j++) {
            ret = fscanf(fh, "%d ", &cn[j]);
            assert(ret == 1);
        }
        ret = fscanf(fh, "%d\n", &cn[stride-1]);
    }

    // CELLS
    fgets(buf, 256, fh);
    fgets(buf, 256, fh);
    next = strtok(buf, "\n");

    bad_ptr = NULL;
    nc = strtod(next, &bad_ptr);
    assert(*bad_ptr == '\0');
    printf("Cells: %d\n", nc);

    C.resize(nc);

    for (E_Int i = 0; i < nc; i++) {
        E_Int stride;
        ret = fscanf(fh, "%d ", &stride);
        assert(ret == 1);
        auto &cn = C[i];
        cn.resize(stride);
        for (E_Int j = 0; j < stride-1; j++) {
            ret = fscanf(fh, "%d ", &cn[j]);
            assert(ret == 1);
        }
        ret = fscanf(fh, "%d\n", &cn[stride-1]);
    }

    fclose(fh);

    make_skin();

    make_bbox();

    hash_skin();

    make_point_faces();

    init_adaptation_data();

    srand(time(NULL));
}

void Mesh::make_bbox()
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

    dmin = std::min(xmin, std::min(ymin, zmin));
    dmax = std::max(xmax, std::max(ymax, zmax));
}

void Mesh::make_skin()
{
    skin.clear();

    std::vector<E_Int> count(nf, 0);
    
    for (const auto &cn : C) {
        for (E_Int face : cn)
            count[face]++;
    }

    for (E_Int i = 0; i < nf; i++) {
        E_Int c = count[i];
        assert(c == 1 || c == 2);
        if (c == 1) skin.push_back(i);
    }
}

void Mesh::hash_skin()
{
    // Throw out the z-coordinate and hash the AABB of the skin faces to a
    // 2D array.
    //std::vector<E_Int, std::vector<E_Int>> box_to_faces;

    NBIN = 100;
    DX = (xmax - xmin);
    DY = (ymax - ymin);
    DZ = (zmax - zmin);

    std::vector<E_Int> pbins(np, -1);

    // Hash the points
    // TODO: hash only the skin points
    for (E_Int i = 0; i < np; i++) {
        E_Float px = X[i];
        E_Float py = Y[i];
        E_Int I = floor(px / DX * NBIN);
        E_Int J = floor(py / DY * NBIN);
        E_Int bin = I + NBIN * J;
        pbins[i] = bin;
    }

    // Hash the faces
    // A face belongs to as many buckets as its points
    bin_faces.clear();

    for (E_Int face : skin) {
        const auto &cn = F[face];
        E_Float xmin, xmax, ymin, ymax;
        xmin = ymin = std::numeric_limits<E_Float>::max();
        xmax = ymax = std::numeric_limits<E_Float>::min();
        for (E_Int p : cn) {
            xmin = std::min(X[p], xmin);
            xmax = std::max(X[p], xmax);
            ymin = std::min(Y[p], ymin);
            ymax = std::max(Y[p], ymax);
        }

        E_Int Imin = floor(xmin / DX * NBIN);
        E_Int Imax = floor(xmax / DX * NBIN);
        E_Int Jmin = floor(ymin / DY * NBIN);
        E_Int Jmax = floor(ymax / DY * NBIN);

        for (E_Int J = Jmin; J < Jmax; J++) {
            for (E_Int I = Imin; I < Imax; I++) {
                E_Int bin = I + NBIN * J;
                bin_faces[bin].insert(face);
            }
        }
    }

    /*
    for (auto &bdata : bin_faces) {
        printf("%d -> ", bdata.first);
        for (E_Int face : bdata.second)
            printf("%d ", face);
        puts("");
    }
    */

}

bool Mesh::is_point_inside(E_Float px, E_Float py, E_Float pz)
{
    // point must be in bounding box
    if (!(xmin <= px && px <= xmax &&
          ymin <= py && py <= ymax &&
          zmin <= pz && pz <= zmax))
        return false;
    
    return true;

    // Count the hits
    E_Int hits = 0;

    TriangleIntersection TI;
    E_Int a, b, c, hit;

    // Choose a random ray direction
    E_Float dx = 0;
    E_Float dy = 0;
    E_Float dz = 1;

    for (E_Int face : skin) {
        const auto &cn = F[face];

        // First triangle
        a = cn[0]; b = cn[1]; c = cn[2];

        hit = MollerTrumbore(px, py, pz, dx, dy, dz, X[a], Y[a], Z[a],
            X[b], Y[b], Z[b], X[c], Y[c], Z[c], TI);
        
        if (hit) {
            hits++;
            continue;
        }

        if (face_is_tri(face))
            continue;
        
        // Second triangle

        a = cn[0]; b = cn[2]; c = cn[3];

        hits += MollerTrumbore(px, py, pz, dx, dy, dz, X[a], Y[a], Z[a],
            X[b], Y[b], Z[b], X[c], Y[c], Z[c], TI);
    }

    // point is inside if number of intersections is odd
    return hits % 2 == 1;
}

void Mesh::write_ngon(const char *fname)
{
    FILE *fh = fopen(fname, "w");
    assert(fh);

    fprintf(fh, "POINTS\n");
    fprintf(fh, "%d\n", np);
    for (E_Int i = 0; i < np; i++) {
        fprintf(fh, "%f %f %f\n", X[i], Y[i], Z[i]);
    }

    fprintf(fh, "INDPG\n");
    fprintf(fh, "%d\n", nf+1);
    E_Int sizeNGon = 0;
    fprintf(fh, "%d ", sizeNGon);
    for (E_Int i = 0; i < nf; i++) {
        sizeNGon += F[i].size();
        fprintf(fh, "%d ", sizeNGon);
    }
    fprintf(fh, "\n");

    fprintf(fh, "NGON\n");
    fprintf(fh, "%d\n", sizeNGon);
    for (E_Int i = 0; i < nf; i++) {
        for (E_Int p : F[i])
            fprintf(fh, "%d ", p);
    }
    fprintf(fh, "\n");

    fprintf(fh, "INDPH\n");
    fprintf(fh, "%d\n", nc+1);
    E_Int sizeNFace = 0;
    fprintf(fh, "%d ", sizeNFace);
    for (E_Int i = 0; i < nc; i++) {
        sizeNFace += C[i].size();
        fprintf(fh, "%d ", sizeNFace);
    }
    fprintf(fh, "\n");

    fprintf(fh, "NFace\n");
    fprintf(fh, "%d\n", sizeNFace);
    for (E_Int i = 0; i < nc; i++) {
        for (E_Int p : C[i])
            fprintf(fh, "%d ", p);
    }
    fprintf(fh, "\n");

    fclose(fh); 
}


bool Mesh::face_contains_sface(E_Int face, E_Int sface, const Mesh &S) const
{
    // face containes mface iff it contains all its points
    assert(S.face_is_active(sface));
    const auto &cn = S.F[sface];
    for (E_Int p : cn) {
        if (face_contains_point(face, S.X[p], S.Y[p]) == -1) return false;
    }
    return true;
}

E_Int Mesh::face_contains_point(E_Int face, E_Float x, E_Float y) const
{
    const auto &cn = F[face];

    E_Int hit, a, b, c;

    if (face_is_quad(face)) {

        // First triangle 

        a = cn[0], b = cn[1], c = cn[2];
        
        hit = Triangle::ispointInside(x, y, X[a], Y[a], X[b], Y[b], X[c], Y[c]);
        
        if (hit) return 0;

        // Second triangle
        a = cn[0], b = cn[2], c = cn[3];
        
        hit = Triangle::ispointInside(x, y, X[a], Y[a], X[b], Y[b], X[c], Y[c]);
        
        if (hit) return 1;
    } else {

        assert(face_is_tri(face));

        a = cn[0], b = cn[1], c = cn[2];
        
        hit = Triangle::ispointInside(x, y, X[a], Y[a], X[b], Y[b], X[c], Y[c]);
        
        if (hit) return 0;
    }

    return -1;
}

UEdge::UEdge(E_Int P, E_Int Q)
{
    p = std::min(P, Q);
    q = std::max(P, Q);
}

bool UEdge::operator<(const UEdge &E) const
{
    return (p < E.p) || (p == E.p && q < E.q);
}

Mesh Mesh::extract_conformized()
{
    // Keep all the points
    std::vector<E_Float> new_X(X), new_Y(Y), new_Z(Z);

    // Conformize the faces

    std::vector<std::vector<E_Int>> new_F(factive.size());

    E_Int new_nf = 0;
    
    std::map<E_Int, E_Int> new_fids;

    for (E_Int face : factive) {
        new_fids[face] = new_nf;

        const auto &pn = F[face];

        auto &new_face = new_F[new_nf];

        for (size_t j = 0; j < pn.size(); j++) {
            E_Int p = pn[j];
            E_Int q = pn[(j+1)%pn.size()];

            UEdge e(p, q);

            auto it = ecenter.find(e);

            if (it != ecenter.end()) {
                new_face.push_back(p);
                new_face.push_back(it->second);
            } else {
                new_face.push_back(p);
            }
        }

        new_nf++;
    }

    // Update cell connectivity

    std::vector<std::vector<E_Int>> new_C(C.size());

    for (E_Int i = 0; i < nc; i++) {
        const auto &pf = C[i];

        auto &new_cell = new_C[i];

        for (E_Int face : pf) {

            if (face_is_active(face)) {
                new_cell.push_back(new_fids[face]);
            } else {
                std::vector<E_Int> fleaves;
                get_fleaves(face, fleaves);

                for (E_Int fleaf : fleaves)
                    new_cell.push_back(new_fids[fleaf]);
            }
        }
    }

    Mesh new_M;
    new_M.np = np;
    new_M.X = X;
    new_M.Y = Y;
    new_M.Z = Z;
    new_M.nf = new_nf;
    new_M.F = new_F;
    new_M.nc = nc;
    new_M.C = new_C;

    for (E_Int face : patch) {
        new_M.patch.insert(new_fids[face]);
    }

    return new_M;
}

void Mesh::get_fleaves(E_Int face, std::vector<E_Int> &fleaves)
{
    if (face_is_active(face)) {
        fleaves.push_back(face);
        return;
    }

    for (E_Int child : fchildren.at(face)) get_fleaves(child, fleaves);
}

PyObject *Mesh::export_karray()
{
    E_Int sizeNGon = 0, sizeNFace = 0;

    for (const auto &pn : F) sizeNGon += (E_Int)pn.size();
    for (const auto &pf : C) sizeNFace += (E_Int)pf.size();

    const char *varString = "CoordinateX,CoordinateY,CoordinateZ";

    PyObject *ret = K_ARRAY::buildArray3(3, varString, np, nc, nf, "NGON",
        sizeNGon, sizeNFace, 3, false, 3);
    
    K_FLD::FldArrayF *f;
    K_FLD::FldArrayI *cn;
    K_ARRAY::getFromArray3(ret, f, cn);

    E_Float *px = f->begin(1);
    for (E_Int i = 0; i < np; i++) px[i] = X[i];
    E_Float *py = f->begin(2);
    for (E_Int i = 0; i < np; i++) py[i] = Y[i];
    E_Float *pz = f->begin(3);
    for (E_Int i = 0; i < np; i++) pz[i] = Z[i];

    E_Int *indPG = cn->getIndPG();
    E_Int *ngon = cn->getNGon();
    E_Int *indPH = cn->getIndPH();
    E_Int *nface = cn->getNFace();

    indPG[0] = indPH[0] = 0;
    for (E_Int i = 0; i < nf; i++) indPG[i+1] = indPG[i] + (E_Int)F[i].size();
    for (E_Int i = 0; i < nc; i++) indPH[i+1] = indPH[i] + (E_Int)C[i].size();

    assert(indPG[nf] == sizeNGon);
    assert(indPH[nc] == sizeNFace);

    E_Int *ptr = ngon;

    for (E_Int i = 0; i < nf; i++) {
        const auto &pn = F[i];
        for (E_Int p : pn) *ptr++ = p+1;
    }

    ptr = nface;

    for (E_Int i = 0; i < nc; i++) {
        const auto &pf = C[i];
        for (E_Int f : pf) *ptr++ = f+1;
    }

    return ret;
}
