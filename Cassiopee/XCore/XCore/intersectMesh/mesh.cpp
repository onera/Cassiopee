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
    Int p, q;

    DEdge(Int P, Int Q) : p(P), q(Q) {}

    bool operator<(const DEdge &f) const
    {
        Int ep = std::min(p, q);
        Int eq = std::max(p, q);
        Int fp = std::min(f.p, f.q);
        Int fq = std::max(f.p, f.q);
        return (ep < fp) || (ep == fp && eq < fq);
    }
};

void IMesh::make_edges()
{
    std::map<DEdge, Int> edges;

    ne = 0;

    F2E.resize(nf);

    for (Int i = 0; i < nf; i++) {
        const auto &pn = F[i];
        
        F2E[i].resize(pn.size());

        for (size_t j = 0; j < pn.size(); j++) {
            Int p = pn[j];
            Int q = pn[(j+1)%pn.size()];
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


std::vector<pointFace> IMesh::locate(Int p, Float px, Float py, Float pz,
    const std::set<Int> &patch) const
{
    Int a, b, c;
    Int hit;
    std::vector<pointFace> hits;
    
    for (Int face : patch) {
        assert(face_is_active(face));

        const auto &cn = F[face];

        // First triangle

        a = cn[0], b = cn[1], c = cn[2];

        hit = Triangle::isPointInside(px, py, pz, X[a], Y[a], Z[a],
            X[b], Y[b], Z[b], X[c], Y[c], Z[c]);

        if (hit) {
            hits.push_back(pointFace(face, 0));
            continue;
        }
        
        if (face_is_tri(face)) continue;

        // Second triangle
        a = cn[0], b = cn[2], c = cn[3];

        hit = Triangle::isPointInside(px, py, pz, X[a], Y[a], Z[a],
            X[b], Y[b], Z[b], X[c], Y[c], Z[c]);

        if (hit) {
            hits.push_back(pointFace(face, 1));
        }
    }

    return hits;
}

void IMesh::init_adaptation_data()
{
    flevel.resize(nf, 0);

    for (Int i = 0; i < nf; i++) factive.insert(i);
}

bool IMesh::faces_are_dups(Int mface, Int sface, const IMesh &S)
{
    const auto &pnm = F[mface];
    const auto &pns = S.F[sface];

    assert(face_is_quad(mface) || face_is_tri(mface));
    assert(S.face_is_quad(sface) || S.face_is_tri(sface));

    if (pnm.size() != pns.size()) return false;

    Int mfound[4] = { 0, 0, 0, 0 };

    for (size_t i = 0; i < pnm.size(); i++) {
        Int pm = pnm[i];
        for (size_t j = 0; j < pns.size(); j++) {
            Int ps = pns[j];
            if (cmp_points(X[pm], Y[pm], Z[pm], S.X[ps], S.Y[ps], S.Z[ps]) == 0) {
                assert(mfound[i] == 0);
                mfound[i] = 1;
                break;
            }
        }

        if (mfound[i] == 0) return false;
    }

    return true;
}

IMesh::IMesh()
{}

IMesh::IMesh(K_FLD::FldArrayI &cn, Float *x, Float *y, Float *z, Int npts)
{
    np = npts;
    ne = 0;
    nf = cn.getNFaces();
    nc = cn.getNElts();

    X.resize(np);
    Y.resize(np);
    Z.resize(np);
    for (Int i = 0; i < np; i++) {
        X[i] = x[i];
        Y[i] = y[i];
        Z[i] = z[i];
    }

    F.reserve(nf);
    Int *indPG = cn.getIndPG();
    Int *ngon = cn.getNGon();
    for (Int i = 0; i < nf; i++) {
        Int np = -1;
        Int *pn = cn.getFace(i, np, ngon, indPG);
        std::vector<Int> points(np);
        for (Int j = 0; j < np; j++)
            points[j] = pn[j] - 1;
        F.push_back(points);
    }

    C.reserve(nc);
    Int *indPH = cn.getIndPH();
    Int *nface = cn.getNFace();
    for (Int i = 0; i < nc; i++) {
        Int nf = -1;
        Int *pf = cn.getElt(i, nf, nface, indPH);
        std::vector<Int> faces(nf);
        for (Int j = 0; j < nf; j++)
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

void IMesh::make_point_faces()
{
    P2F.clear();
    P2F.resize(np);
    
    for (Int face : factive) {
        const auto &pn = F[face];
        for (Int p : pn) P2F[p].push_back(face);
    }
}

IMesh::IMesh(const char *fname)
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
    printf("points: " SF_D_ "\n", np);

    X.resize(np);
    Y.resize(np);
    Z.resize(np);

    int ret;

    for (Int i = 0; i < np; i++) {
        ret = fscanf(fh, "%lf %lf %lf\n", &X[i], &Y[i], &Z[i]);
        if (ret != 3) abort();
        //assert(ret == 3);
    }

    // FACES
    fgets(buf, 256, fh);
    fgets(buf, 256, fh);
    next = strtok(buf, "\n");

    bad_ptr = NULL;
    nf = strtod(next, &bad_ptr);
    assert(*bad_ptr == '\0');
    printf("Faces: " SF_D_ "\n", nf);

    F.resize(nf);

    for (Int i = 0; i < nf; i++) {
        Int stride;
        ret = fscanf(fh, SF_D_ " ", &stride);
        if (ret != 1) abort();
        auto &cn = F[i];
        cn.resize(stride);
        for (Int j = 0; j < stride-1; j++) {
            ret = fscanf(fh, SF_D_ " ", &cn[j]);
            if (ret != 1) abort();
        }
        ret = fscanf(fh, SF_D_ "\n", &cn[stride-1]);
        if (ret != 1) abort();
    }

    // CELLS
    fgets(buf, 256, fh);
    fgets(buf, 256, fh);
    next = strtok(buf, "\n");

    bad_ptr = NULL;
    nc = strtod(next, &bad_ptr);
    assert(*bad_ptr == '\0');
    printf("Cells: " SF_D_ "\n", nc);

    C.resize(nc);

    for (Int i = 0; i < nc; i++) {
        Int stride;
        ret = fscanf(fh, SF_D_ " ", &stride);
        if (ret != 1) abort();
        auto &cn = C[i];
        cn.resize(stride);
        for (Int j = 0; j < stride-1; j++) {
            ret = fscanf(fh, SF_D_ " ", &cn[j]);
            if (ret != 1) abort();
        }
        ret = fscanf(fh, SF_D_ "\n", &cn[stride-1]);
        if (ret != 1) abort();
    }

    fclose(fh);

    make_skin();

    make_bbox();

    hash_skin();

    make_point_faces();

    init_adaptation_data();

    srand(time(NULL));
}

void IMesh::make_bbox()
{
    xmin = ymin = zmin = std::numeric_limits<Float>::max();
    xmax = ymax = zmax = std::numeric_limits<Float>::min();

    for (Int i = 0; i < np; i++) {
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

void IMesh::make_skin()
{
    skin.clear();

    std::vector<Int> count(nf, 0);
    
    for (const auto &cn : C) {
        for (Int face : cn)
            count[face]++;
    }

    for (Int i = 0; i < nf; i++) {
        Int c = count[i];
        assert(c == 1 || c == 2);
        if (c == 1) skin.push_back(i);
    }
}

void IMesh::hash_skin()
{
    // Throw out the z-coordinate and hash the AABB of the skin faces to a
    // 2D array.
    //std::vector<Int, std::vector<Int>> box_to_faces;

    NBIN = 100;
    DX = (xmax - xmin);
    DY = (ymax - ymin);
    DZ = (zmax - zmin);

    std::vector<Int> pbins(np, -1);

    // Hash the points
    // TODO: hash only the skin points
    for (Int i = 0; i < np; i++) {
        Float px = X[i];
        Float py = Y[i];
        Int I = floor(px / DX * NBIN);
        Int J = floor(py / DY * NBIN);
        Int bin = I + NBIN * J;
        pbins[i] = bin;
    }

    // Hash the faces
    // A face belongs to as many buckets as its points
    bin_faces.clear();

    for (Int face : skin) {
        const auto &cn = F[face];
        Float xmin, xmax, ymin, ymax;
        xmin = ymin = std::numeric_limits<Float>::max();
        xmax = ymax = std::numeric_limits<Float>::min();
        for (Int p : cn) {
            xmin = std::min(X[p], xmin);
            xmax = std::max(X[p], xmax);
            ymin = std::min(Y[p], ymin);
            ymax = std::max(Y[p], ymax);
        }

        Int Imin = floor(xmin / DX * NBIN);
        Int Imax = floor(xmax / DX * NBIN);
        Int Jmin = floor(ymin / DY * NBIN);
        Int Jmax = floor(ymax / DY * NBIN);

        for (Int J = Jmin; J < Jmax; J++) {
            for (Int I = Imin; I < Imax; I++) {
                Int bin = I + NBIN * J;
                bin_faces[bin].insert(face);
            }
        }
    }
}

bool IMesh::is_point_inside(Float px, Float py, Float pz)
{
    // point must be in bounding box
    if (!(xmin <= px && px <= xmax &&
          ymin <= py && py <= ymax &&
          zmin <= pz && pz <= zmax))
        return false;
    
    // Count the hits
    Int hits = 0;

    TriangleIntersection TI;
    Int a, b, c, hit;

    // Choose a random ray direction
    Float dx = 0.2;
    Float dy = -0.5;
    Float dz = 0.4;

    for (Int face : skin) {
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

void IMesh::write_ngon(const char *fname)
{
    FILE *fh = fopen(fname, "w");
    assert(fh);

    fprintf(fh, "POINTS\n");
    fprintf(fh, SF_D_ "\n", np);
    for (Int i = 0; i < np; i++) {
        fprintf(fh, "%f %f %f\n", X[i], Y[i], Z[i]);
    }

    fprintf(fh, "INDPG\n");
    fprintf(fh, SF_D_ "\n", nf+1);
    Int sizeNGon = 0;
    fprintf(fh, SF_D_ " ", sizeNGon);
    for (Int i = 0; i < nf; i++) {
        sizeNGon += F[i].size();
        fprintf(fh, SF_D_ " ", sizeNGon);
    }
    fprintf(fh, "\n");

    fprintf(fh, "NGON\n");
    fprintf(fh, SF_D_ "\n", sizeNGon);
    for (Int i = 0; i < nf; i++) {
        for (Int p : F[i])
            fprintf(fh, SF_D_ " ", p);
    }
    fprintf(fh, "\n");

    fprintf(fh, "INDPH\n");
    fprintf(fh, SF_D_ "\n", nc+1);
    Int sizeNFace = 0;
    fprintf(fh, SF_D_ " ", sizeNFace);
    for (Int i = 0; i < nc; i++) {
        sizeNFace += C[i].size();
        fprintf(fh, SF_D_ " ", sizeNFace);
    }
    fprintf(fh, "\n");

    fprintf(fh, "NFace\n");
    fprintf(fh, SF_D_ "\n", sizeNFace);
    for (Int i = 0; i < nc; i++) {
        for (Int p : C[i])
            fprintf(fh, SF_D_ " ", p);
    }
    fprintf(fh, "\n");

    fclose(fh); 
}


bool IMesh::face_contains_sface(Int face, Int sface, const IMesh &S) const
{
    // face containes mface iff it contains all its points
    assert(S.face_is_active(sface));
    const auto &cn = S.F[sface];
    for (Int p : cn) {
        if (face_contains_point(face, S.X[p], S.Y[p], S.Z[p]) == -1) 
            return false;
    }
    return true;
}

Int IMesh::face_contains_point(Int face, Float x, Float y, Float z) const
{
    const auto &cn = F[face];

    Int hit, a, b, c;

    if (face_is_quad(face)) {

        // First triangle 

        a = cn[0], b = cn[1], c = cn[2];
        
        hit = Triangle::isPointInside(x, y, z, X[a], Y[a], Z[a],
            X[b], Y[b], Z[b], X[c], Y[c], Z[c]);
        
        if (hit) return 0;

        // Second triangle
        a = cn[0], b = cn[2], c = cn[3];
        
        hit = Triangle::isPointInside(x, y, z, X[a], Y[a], Z[a],
            X[b], Y[b], Z[b], X[c], Y[c], Z[c]);
        
        if (hit) return 1;
    } else {

        assert(face_is_tri(face));

        a = cn[0], b = cn[1], c = cn[2];
        
        hit = Triangle::isPointInside(x, y, z, X[a], Y[a], Z[a],
            X[b], Y[b], Z[b], X[c], Y[c], Z[c]);
        
        if (hit) return 0;
    }

    return -1;
}

UEdge::UEdge(Int P, Int Q)
{
    p = std::min(P, Q);
    q = std::max(P, Q);
}

bool UEdge::operator<(const UEdge &E) const
{
    return (p < E.p) || (p == E.p && q < E.q);
}

IMesh IMesh::extract_conformized()
{
    // Keep all the points
    std::vector<Float> new_X(X), new_Y(Y), new_Z(Z);

    // Conformize the faces

    std::vector<std::vector<Int>> new_F(factive.size());

    Int new_nf = 0;
    
    std::map<Int, Int> new_fids;

    for (Int face : factive) {
        new_fids[face] = new_nf;

        const auto &pn = F[face];

        auto &new_face = new_F[new_nf];

        for (size_t j = 0; j < pn.size(); j++) {
            Int p = pn[j];
            Int q = pn[(j+1)%pn.size()];

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

    std::vector<std::vector<Int>> new_C(C.size());

    for (Int i = 0; i < nc; i++) {
        const auto &pf = C[i];

        auto &new_cell = new_C[i];

        for (Int face : pf) {

            if (face_is_active(face)) {
                new_cell.push_back(new_fids[face]);
            } else {
                std::vector<Int> fleaves;
                get_fleaves(face, fleaves);

                for (Int fleaf : fleaves)
                    new_cell.push_back(new_fids[fleaf]);
            }
        }
    }

    IMesh new_M;
    new_M.np = np;
    new_M.X = X;
    new_M.Y = Y;
    new_M.Z = Z;
    new_M.nf = new_nf;
    new_M.F = new_F;
    new_M.nc = nc;
    new_M.C = new_C;

    for (Int face : patch) {
        new_M.patch.insert(new_fids[face]);
    }

    return new_M;
}

void IMesh::get_fleaves(Int face, std::vector<Int> &fleaves)
{
    if (face_is_active(face)) {
        fleaves.push_back(face);
        return;
    }

    for (Int child : fchildren.at(face)) get_fleaves(child, fleaves);
}

PyObject *IMesh::export_karray()
{
    Int sizeNGon = 0, sizeNFace = 0;

    for (const auto &pn : F) sizeNGon += (Int)pn.size();
    for (const auto &pf : C) sizeNFace += (Int)pf.size();

    const char *varString = "CoordinateX,CoordinateY,CoordinateZ";

    PyObject *ret = K_ARRAY::buildArray3(3, varString, np, nc, nf, "NGON",
        sizeNGon, sizeNFace, 3, false, 3);
    
    K_FLD::FldArrayF *f;
    K_FLD::FldArrayI *cn;
    K_ARRAY::getFromArray3(ret, f, cn);

    Float *px = f->begin(1);
    for (Int i = 0; i < np; i++) px[i] = X[i];
    Float *py = f->begin(2);
    for (Int i = 0; i < np; i++) py[i] = Y[i];
    Float *pz = f->begin(3);
    for (Int i = 0; i < np; i++) pz[i] = Z[i];

    Int *indPG = cn->getIndPG();
    Int *ngon = cn->getNGon();
    Int *indPH = cn->getIndPH();
    Int *nface = cn->getNFace();

    indPG[0] = indPH[0] = 0;
    for (Int i = 0; i < nf; i++) indPG[i+1] = indPG[i] + (Int)F[i].size();
    for (Int i = 0; i < nc; i++) indPH[i+1] = indPH[i] + (Int)C[i].size();

    assert(indPG[nf] == sizeNGon);
    assert(indPH[nc] == sizeNFace);

    Int *ptr = ngon;

    for (Int i = 0; i < nf; i++) {
        const auto &pn = F[i];
        for (Int p : pn) *ptr++ = p+1;
    }

    ptr = nface;

    for (Int i = 0; i < nc; i++) {
        const auto &pf = C[i];
        for (Int f : pf) *ptr++ = f+1;
    }

    delete f;
    delete cn;


    return ret;
}
