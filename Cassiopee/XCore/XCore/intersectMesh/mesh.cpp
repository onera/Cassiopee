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
#include "io.h"

void IMesh::triangulate_face_set()
{
    E_Int NF = nf;

    E_Int face_incr = (E_Int)faces_to_tri.size();

    F.resize(NF + face_incr);

    std::vector<E_Int> owner(NF, -1), neigh(NF, -1);

    for (E_Int i = 0; i < nc; i++) {
        const auto &pf = C[i];
        for (E_Int fid : pf) {
            if (owner[fid] == -1) owner[fid] = i;
            else neigh[fid] = i;
        }
    }

    assert(skin.size() > 0);

    for (E_Int fid : faces_to_tri) {

        auto &pn = F[fid];

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

        NF++;
    }

    assert(NF == nf + face_incr);

    nf = NF;

    faces_to_tri.clear();

    for (E_Int cid = 0; cid < nc; cid++) {
        assert(C[cid].size() == 6 || C[cid].size() == 7);
    }
}

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

void IMesh::make_edges()
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

void IMesh::init_adaptation_data()
{
    flevel.resize(nf, 0);

    for (E_Int i = 0; i < nf; i++) factive.insert(i);
}

bool IMesh::faces_are_dups(E_Int mface, E_Int sface, const IMesh &S)
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

IMesh::IMesh(K_FLD::FldArrayI &cn, E_Float *x, E_Float *y, E_Float *z, E_Int npts)
{
    NX = 100;
    NY = 100;
    NZ = 100;
    NXY = NX * NY;
    NXYZ = NXY * NZ;

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

    //make_skin();

    //make_bbox();

    //hash_skin();

    //make_point_faces();

    //init_adaptation_data();
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

E_Int IMesh::RayFaceIntersect(E_Float px, E_Float py, E_Float pz, E_Float dx,
    E_Float dy, E_Float dz, E_Int fid, TriangleIntersection &TI) const
{
    const auto &pn = F[fid];

    // TODO(Imad): hexa mesh for now
    assert(pn.size() == 4);

    E_Int a = pn[0], b = pn[1], c = pn[2];

    E_Int hit;

    TI.face = fid;

    TI.tri = 0;

    hit = MollerTrumbore(px, py, pz, dx, dy, dz, X[a], Y[a], Z[a], X[b], Y[b],
        Z[b], X[c], Y[c], Z[c], TI);
    
    if (hit) {
        if      (Sign(TI.u) == 0 && Sign(TI.v) == 0) TI.vid = 0;
        else if (Sign(TI.v) == 0 && Sign(TI.w) == 0) TI.vid = 1;
        else if (Sign(TI.w) == 0 && Sign(TI.u) == 0) TI.vid = 2;
        else if (Sign(TI.v) == 0) TI.eid = 0;
        else if (Sign(TI.w) == 0) TI.eid = 1;

        return 1;
    }

    E_Int d = pn[3];

    TI.tri = 1;

    hit = MollerTrumbore(px, py, pz, dx, dy, dz, X[c], Y[c], Z[c], X[d], Y[d],
        Z[d], X[a], Y[a], Z[a], TI);

    if (hit) {
        if      (Sign(TI.u) == 0 && Sign(TI.v) == 0) TI.vid = 2;
        else if (Sign(TI.v) == 0 && Sign(TI.w) == 0) TI.vid = 3;
        else if (Sign(TI.w) == 0 && Sign(TI.u) == 0) TI.vid = 0;
        else if (Sign(TI.v) == 0) TI.eid = 2;
        else if (Sign(TI.w) == 0) TI.eid = 3;

        return 1;
    }

    return 0;
}

E_Int IMesh::project_point(E_Float ox, E_Float oy, E_Float oz, E_Float dx,
    E_Float dy, E_Float dz, TriangleIntersection &TI, E_Int II)
{
    // Calculate entry point

    E_Float tentry = 0;
    E_Float texit = 0;

    E_Float px = ox;
    E_Float py = oy;
    E_Float pz = oz;

    if (ox < xmin || ox > xmax ||
        oy < ymin || oy > ymax ||
        oz < zmin || oz > zmax) {

        // Ray origin outside mesh bounding box, compute the intersection

        E_Float tminx, tmaxx, tminy, tmaxy, tminz, tmaxz;

        tminx = (xmin - ox) / dx;
        tmaxx = (xmax - ox) / dx;
        if (dx < 0) std::swap(tminx, tmaxx);

        tminy = (ymin - oy) / dy;
        tmaxy = (ymax - oy) / dy;
        if (dy < 0) std::swap(tminy, tmaxy);

        tminz = (zmin - oz) / dz;
        tmaxz = (zmax - oz) / dz;
        if (dz < 0) std::swap(tminz, tmaxz);

        tentry = std::max(tminx, std::max(tminy, tminz));
        texit = std::min(tmaxx, std::min(tmaxy, tmaxz));

        // Check for intersection

        if (tentry > texit || texit < 0) return 0;

        px = ox + tentry * dx;
        py = oy + tentry * dy;
        pz = oz + tentry * dz;

        // Make sure the entry lies within bbox

        if (Sign(px - xmin) == 0) px = xmin;
        if (Sign(py - ymin) == 0) py = ymin;
        if (Sign(pz - zmin) == 0) pz = zmin;
        if (Sign(px - xmax) == 0) px = xmax;
        if (Sign(py - ymax) == 0) py = ymax;
        if (Sign(pz - zmax) == 0) pz = zmax;
    }

    //point_write("entry", px, py, pz);

    ox = px, oy = py, oz = pz;

    E_Int voxel_x = floor((ox - xmin) / HX);
    E_Int voxel_y = floor((oy - ymin) / HY);
    E_Int voxel_z = floor((oz - zmin) / HZ);

    assert(voxel_x >= 0 && voxel_x < NX);
    assert(voxel_y >= 0 && voxel_y < NY);
    assert(voxel_z >= 0 && voxel_z < NZ);

    // Steps

    E_Int step_x, step_y, step_z;
    E_Float tDeltaX, tDeltaY, tDeltaZ;
    E_Float tmaxX, tmaxY, tmaxZ;

    if (dx > 0) {
        step_x = 1;
        tDeltaX = HX/dx;
        tmaxX = ((voxel_x+1)*HX + xmin - ox) / dx;
    } else if (dx < 0) {
        step_x = -1;
        tDeltaX = -HX/dx;
        tmaxX = (voxel_x*HX + xmin - ox) / dx;
    } else {
        step_x = 0;
        tDeltaX = EFLOATMAX;
        tmaxX = EFLOATMAX;
    }

    if (dy > 0) {
        step_y = 1;
        tDeltaY = HY/dy;
        tmaxY = ((voxel_y+1)*HY + ymin - oy) / dy;
    } else if (dy < 0) {
        step_y = -1;
        tDeltaY = -HY/dy;
        tmaxY = (voxel_y*HY + ymin - oy) / dy;
    } else {
        step_y = 0;
        tDeltaY = EFLOATMAX;
        tmaxY = EFLOATMAX;
    }
    
    if (dz > 0) {
        step_z = 1;
        tDeltaZ = HZ/dz;
        tmaxZ = ((voxel_z+1)*HZ + zmin - oz) / dz;
    } else if (dz < 0) {
        step_z = -1;
        tDeltaZ = -HZ/dz;
        tmaxZ = (voxel_z*HZ + zmin - oz) / dz;
    } else {
        step_z = 0;
        tDeltaZ = EFLOATMAX;
        tmaxZ = EFLOATMAX;
    }

    assert(tDeltaX >= 0);
    assert(tDeltaY >= 0);
    assert(tDeltaZ >= 0);
    assert(tmaxX >= 0);
    assert(tmaxY >= 0);
    assert(tmaxZ >= 0);

    //edge_write("bad_edge", ox, oy, oz, ox+100*dx, oy+100*dy, oz+100*dz);

    E_Int current_cell = get_voxel(voxel_x, voxel_y, voxel_z);

    std::set<E_Int> faces_set;

    std::set<E_Int> tested_faces;

    TriangleIntersection hitTI;

    E_Int hit = 0;

    while (current_cell < NXYZ && current_cell >= 0) {

        // Check for intersections in the current cell

        const auto &pf = bin_faces[current_cell];

        for (auto fid : pf) faces_set.insert(fid);

        for (E_Int fid : pf) {
            if (tested_faces.find(fid) != tested_faces.end()) continue;
            tested_faces.insert(fid);
            hit += RayFaceIntersect(ox, oy, oz, dx, dy, dz, fid, hitTI);
            if (hit == 1) {
                TI = hitTI;
                return 1;
            }
        }

        if (tmaxX <= tmaxY && tmaxX <= tmaxZ) {
            tmaxX = tmaxX + tDeltaX;
            voxel_x += step_x;
        } else if (tmaxY <= tmaxX && tmaxY <= tmaxZ) {
            tmaxY = tmaxY + tDeltaY;
            voxel_y += step_y;
        } else if (tmaxZ <= tmaxX && tmaxZ <= tmaxY) {
            tmaxZ = tmaxZ + tDeltaZ;
            voxel_z += step_z;
        } else {
            assert(0);
        }

        current_cell = get_voxel(voxel_x, voxel_y, voxel_z);
    }

    /*
    if (II == 8) {
    puts("WRITING FACES FOR THIS POINT");
    printf("Tested %zu faces\n", faces_set.size());
    std::vector<E_Int> faces_out;
    for (E_Int fid : faces_set) faces_out.push_back(fid);
    write_faces("faces", faces_out);
    }
    */

    return hit > 0;
}

bool IMesh::is_point_inside(E_Float ox, E_Float oy, E_Float oz) const
{
    // point must be in bounding box
    if (xmin > ox || ox > xmax ||
        ymin > oy || oy > ymax ||
        zmin > oz || oz > zmax) {
        return false;
    }

    // Choose a random ray direction
    E_Float dx = 0.0;
    E_Float dy = 0.0;
    E_Float dz = 1.0;

    E_Int NORM = sqrt(dx*dx + dy*dy + dz*dz);
    dx /= NORM, dy /= NORM, dz /= NORM;

    E_Float Ogrid[3] = {ox - xmin, oy - ymin, oz - zmin};
    E_Float Ocell[3] = {Ogrid[0] / HX,
                        Ogrid[1] / HY,
                        Ogrid[2] / HZ};

    E_Int voxel_x = floor(Ocell[0]);
    E_Int voxel_y = floor(Ocell[1]);
    E_Int voxel_z = floor(Ocell[2]);

    assert(voxel_x >= 0 && voxel_x < NX);
    assert(voxel_y >= 0 && voxel_y < NY);
    assert(voxel_z >= 0 && voxel_z < NZ);

    E_Float tx, ty, tz;
    tx = ty = tz = EFLOATMAX;

    E_Float deltaTx, deltaTy, deltaTz;
    deltaTx = deltaTy = deltaTz = EFLOATMAX;

    if (dx > 0) {
        tx = (floor(Ocell[0]) + 1) * HX - Ogrid[0];
        tx /= dx;
        deltaTx = HX/dx;
    } else if (dx < 0) {
        tx = floor(Ocell[0]) * HX - Ogrid[0];
        tx /= dx;
        deltaTy = -HX/dx;
    }

    if (dy > 0) {
        ty = (floor(Ocell[1]) + 1) * HY - Ogrid[1];
        ty /= dy;
        deltaTy = HY/dy;
    } else if (dy < 0) {
        ty = floor(Ocell[1]) * HY - Ogrid[1];
        ty /= dy;
        deltaTy = -HY/dy;
    }

    if (dz > 0) {
        tz = (floor(Ocell[2]) + 1) * HZ - Ogrid[2];
        tz /= dz;
        deltaTz = HZ/dz;
    } else if (dz < 0) {
        tz = floor(Ocell[2]) * HZ - Ogrid[2];
        tz /= dz;
        deltaTz = -HZ/dz;
    }

    //E_Float t = 0;

    E_Int current_cell = get_voxel(voxel_x, voxel_y, voxel_z);
    assert(current_cell >= 0);
    assert(current_cell < NXYZ);

    E_Int hits = 0;

    TriangleIntersection TI;

    std::set<E_Int> tested_faces;

    std::set<UEdge> hit_edges;
    std::set<E_Int> hit_vertices;

    // Each face must be tested once
    // Each hit edge must be counted only once
    // Each hit vertex must be counted only once

    while (current_cell < NXYZ && current_cell >= 0) {

        // Check for intersections in the current cell

        //if (it != bin_faces.end()) {

            const auto &pf = bin_faces[current_cell];

            for (E_Int fid : pf) {
                if (tested_faces.find(fid) != tested_faces.end()) continue;
                tested_faces.insert(fid);
                E_Int hit = RayFaceIntersect(ox, oy, oz, dx, dy, dz, fid, TI);
                
                if (hit && TI.t > 0) {

                    const auto &pn = F[fid];

                    // Hit an edge, count it only once
                    if (TI.eid != -1) {
                        assert(TI.vid == -1);

                        E_Int idx = TI.eid;
                        E_Int p = pn[idx];
                        E_Int q = pn[(idx+1)%pn.size()];

                        UEdge E(p, q);

                        if (hit_edges.find(E) == hit_edges.end()) {
                            hit_edges.insert(E);
                            hits++;
                        }
                    }
                    
                    // Hit a vertex, count it only once
                    else if (TI.vid != -1) {
                        E_Int idx = TI.vid;
                        E_Int v = pn[idx];

                        if (hit_vertices.find(v) == hit_vertices.end()) {
                            hit_vertices.insert(v);
                            hits++;
                        }
                    }

                    // Hit the interior of a newly visited face, count it
                    else {
                        assert(TI.vid == -1);
                        assert(TI.eid == -1);

                        hits++;
                    }
                }
            }
        //}

        if (tx < ty && tx && tz) {

            //t = tx;
            tx += deltaTx;
            voxel_x += (dx > 0) ? 1 : -1;

        } else if (ty < tz) {
            assert(ty < tx);

            //t = ty;
            ty += deltaTy;
            voxel_y += (dy > 0) ? 1 : -1;


        } else {
            assert(tz < tx && tz < ty);

            //t = tz;
            tz += deltaTz;
            voxel_z += (dz > 0) ? 1 : -1;

        }

        current_cell = get_voxel(voxel_x, voxel_y, voxel_z);
    }

    // Point is inside if the number of intersections is odd

    /*
    if (hits > 2) {

        E_Float qx = ox + 10000*dx;
        E_Float qy = oy + 10000*dy;
        E_Float qz = oz + 10000*dz;
        edge_write("tested_edge", ox, oy, oz, qx, qy, qz);

        std::vector<E_Int> faces;
        for (auto fid : tested_faces) faces.push_back(fid);
        write_faces("tested_faces", faces);
        printf("oups hits: %d\n", hits);

        assert(0);
    }
    */

    return hits % 2 == 1;
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

void IMesh::write_faces(const char *fname, const std::vector<E_Int> &faces) const
{
    FILE *fh = fopen(fname, "w");
    assert(fh);

    std::map<E_Int, E_Int> new_pids;
    std::map<E_Int, E_Int> old_pids;

    E_Int npts = 0;

    for (E_Int face : faces) {
        const auto &pn = F[face];

        for (E_Int p : pn) {
            if (new_pids.find(p) == new_pids.end()) {
                new_pids[p] = npts;
                old_pids[npts] = p;
                npts++;
            }
        }
    }
    
    fprintf(fh, "POINTS\n");
    fprintf(fh, "%d\n", npts);
    for (E_Int i = 0; i < npts; i++) {
        E_Int opid = old_pids[i];
        fprintf(fh, "%f %f %f\n", X[opid], Y[opid], Z[opid]);
    }

    fclose(fh);
}

void IMesh::write_ngon(const char *fname)
{
    FILE *fh = fopen(fname, "w");
    assert(fh);

    fprintf(fh, "POINTS\n");
    fprintf(fh, SF_D_ "\n", np);
    for (E_Int i = 0; i < np; i++) {
        fprintf(fh, "%f %f %f\n", X[i], Y[i], Z[i]);
    }

    fprintf(fh, "INDPG\n");
    fprintf(fh, SF_D_ "\n", nf+1);
    E_Int sizeNGon = 0;
    fprintf(fh, SF_D_ " ", sizeNGon);
    for (E_Int i = 0; i < nf; i++) {
        sizeNGon += F[i].size();
        fprintf(fh, SF_D_ " ", sizeNGon);
    }
    fprintf(fh, "\n");

    fprintf(fh, "NGON\n");
    fprintf(fh, SF_D_ "\n", sizeNGon);
    for (E_Int i = 0; i < nf; i++) {
        for (E_Int p : F[i])
            fprintf(fh, SF_D_ " ", p);
    }
    fprintf(fh, "\n");

    fprintf(fh, "INDPH\n");
    fprintf(fh, SF_D_ "\n", nc+1);
    E_Int sizeNFace = 0;
    fprintf(fh, SF_D_ " ", sizeNFace);
    for (E_Int i = 0; i < nc; i++) {
        sizeNFace += C[i].size();
        fprintf(fh, SF_D_ " ", sizeNFace);
    }
    fprintf(fh, "\n");

    fprintf(fh, "NFace\n");
    fprintf(fh, SF_D_ "\n", sizeNFace);
    for (E_Int i = 0; i < nc; i++) {
        for (E_Int p : C[i])
            fprintf(fh, SF_D_ " ", p);
    }
    fprintf(fh, "\n");

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