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
#include "triangle.h"
#include "primitives.h"
#include "io.h"

bool is_point_in_2D_polygon(Point2D point, Point2D *polygon, int num_points)
{
    bool inside = false;

    for (int i = 0, j = num_points-1; i < num_points; j = i++) {
        if (((polygon[i].y > point.y) != (polygon[j].y > point.y)) &&
            (point.x < (polygon[j].x - polygon[i].x) * (point.y - polygon[i].y) /
            (polygon[j].y - polygon[i].y) + polygon[i].x)) {
            inside = !inside;
        }
    }
    return inside;
}

Point2D project_to_2D(Point3D point, const E_Float *N)
{
    Point2D projected;

    if (fabs(N[2]) > fabs(N[0]) && fabs(N[2]) > fabs(N[1])) {
        projected.x = point.x;
        projected.y = point.y;
    } else if (fabs(N[1]) > fabs(N[0])) {
        projected.x = point.x;
        projected.y = point.z;
    } else {
        projected.x = point.y;
        projected.y = point.z;
    }

    return projected;
}

#define MAX_PTS 100

bool Smesh::is_point_in_3D_polygon(E_Float x, E_Float y, E_Float z, E_Int fid) const
{
    const auto &pn = Fc[fid];
    const auto *fN = &fnormals[3*fid];
    int a = pn[0];
    E_Float ap[3] = {x-X[a], y-Y[a], z-Z[a]};
    E_Float dp = fabs(K_MATH::dot(ap, fN, 3));
    if (dp > TOL) return false;

    assert(pn.size() <= MAX_PTS);
    Point2D projected_polygon[MAX_PTS];
    for (size_t i = 0; i < pn.size(); i++) {
        Point3D p = {X[pn[i]], Y[pn[i]], Z[pn[i]]};
        projected_polygon[i] = project_to_2D(p, fN); 
    }
    Point2D projected_point = project_to_2D({x, y, z}, fN);


    return is_point_in_2D_polygon(projected_point, projected_polygon, pn.size());
}

bool Smesh::is_point_a_polygon_vertex(E_Float x, E_Float y, E_Float z,
    E_Int fid, PointLoc &ploc, E_Float min_pdist) const
{
    const auto &pn = Fc[fid];

    for (size_t i = 0; i < pn.size(); i++) {
        E_Int p = pn[i];
        E_Float dist = (x-X[p])*(x-X[p]) + (y-Y[p])*(y-Y[p]) + (z-Z[p])*(z-Z[p]);
        dist = sqrt(dist);
        if (dist < min_pdist/10) {
            ploc.v_idx = i;
            ploc.x = X[p];
            ploc.y = Y[p];
            ploc.z = Z[p];
            return true;
        }
    }

    return false;
}

bool Smesh::is_point_on_a_polygon_edge(E_Float x, E_Float y, E_Float z,
    E_Int fid, PointLoc &ploc, E_Float min_pdist) const
{
    const auto &pn = Fc[fid];

    for (size_t i = 0; i < pn.size(); i++) {
        E_Int a = pn[i];
        E_Int b = pn[(i+1)%pn.size()];

        // Check collinearity of ap with ab
        E_Float ap[3] = {x-X[a], y-Y[a], z-Z[a]};
        E_Float ab[3] = {X[b]-X[a], Y[b]-Y[a], Z[b]-Z[a]};
        E_Float C[3];
        K_MATH::cross(ap, ab, C);
        if (fabs(K_MATH::norm(C, 3)) > TOL) continue;
        E_Float dp1 = K_MATH::dot(ap, ab, 3);
        E_Float dp2 = K_MATH::dot(ab, ab, 3);
        if (dp1 >= 0 && dp1 <= dp2) {
            ploc.e_idx = i;
            E_Float t;
            E_Float dx = fabs(X[a]-X[b]);
            E_Float dy = fabs(Y[a]-Y[b]);
            E_Float dz = fabs(Z[a]-Z[b]);
            if (dx > dy && dx > dz) t = (x-X[a])/(X[b]-X[a]);
            else if (dy > dz) t = (y-Y[a])/(Y[b]-Y[a]);
            else t = (z-Z[a])/(Z[b]-Z[a]);
            ploc.x = X[a] + t * (X[b]-X[a]);
            ploc.y = Y[a] + t * (Y[b]-Y[a]);
            ploc.z = Z[a] + t * (Z[b]-Z[a]);
            return true;
        }
    }
    
    return false;
}

std::vector<PointLoc> Smesh::locate2(const Smesh &Sf) const
{
    std::vector<PointLoc> plocs(Sf.np);

    size_t on_vertex = 0, on_edge = 0;

    for (E_Int pid = 0; pid < Sf.np; pid++) {
        E_Float x = Sf.X[pid];
        E_Float y = Sf.Y[pid];
        E_Float z = Sf.Z[pid];
        std::vector<PointLoc> locs;
        ray_intersect_BVH(x, y, z, 0, 1, 0, root_node_idx, locs);
        if (locs.empty()) {
            fprintf(stderr, "Could not locate point %d\n", pid);
            point_write("lost.im", x, y, z);
        }
        assert(locs.size() > 0);
        auto &ploc = plocs[pid];
        E_Float min_abs_t = EFLOATMAX;
        for (const auto &loc : locs) {
            if (fabs(loc.t) < min_abs_t) {
                min_abs_t = fabs(loc.t);
                ploc = loc;
            }
        }
    }
    
    return plocs;
    /*

    //std::vector<Point> oedge, dedge, vpoints;

    for (E_Int pid = 0; pid < Sf.np; pid++) {
        E_Float x = Sf.X[pid];
        E_Float y = Sf.Y[pid];
        E_Float z = Sf.Z[pid];

        E_Int I = floor((x - xmin) / HX);
        E_Int J = floor((y - ymin) / HY);
        E_Int K = floor((z - zmin) / HZ);
        E_Int voxel = get_voxel(I, J, K);

        const auto &pf = bin_faces.at(voxel);

        if (pid == 4941) write_ngon("pf.im", pf);

        bool found = false;

        auto &ploc = plocs[pid];

        for (auto fid : pf) {

            if (pid == 4941) write_face("fid.im", fid);

            found = is_point_in_3D_polygon(x, y, z, fid);

            if (found) {

                //if (pid == 4941) write_face("fid.im", fid);

                ploc.fid = fid;
                
                if (is_point_a_polygon_vertex(x, y, z, fid, ploc, Sf.min_pdist)) {
                    on_vertex++;
                    //vpoints.push_back({x, y, z});
                } else if (is_point_on_a_polygon_edge(x, y, z, fid, ploc, Sf.min_pdist)) {
                    //oedge.push_back({x, y, z});
                    //dedge.push_back({ploc.x, ploc.y, ploc.z});
                    on_edge++;
                }

                break;
            }
        }

        if (!found) {
            fprintf(stderr, "Couldn't locate point %d\n", pid);
            point_write("lost.im", x, y, z);
            write_ngon("bin.im", pf);
        }

        assert(found);
    }

    printf("on vertex: %zu\n", on_vertex);
    printf("on edge: %zu\n", on_edge);

    //point_write("oedge.im", oedge);
    //point_write("dedge.im", dedge);
    //point_write("vpoints.im", vpoints);
    */

    return plocs;
}

std::vector<PointLoc> Smesh::locate(const Smesh &Sf) const
{
    std::vector<PointLoc> ploc(Sf.np);

    for (E_Int pid = 0; pid < Sf.np; pid++) {

        E_Float x = Sf.X[pid];
        E_Float y = Sf.Y[pid];
        E_Float z = Sf.Z[pid];

        E_Int I = floor((x - xmin) / HX);
        E_Int J = floor((y - ymin) / HY);
        E_Int K = floor((z - zmin) / HZ);
        E_Int voxel = get_voxel(I, J, K);

        const auto &pf = bin_faces.at(voxel);

        bool found = false;

        auto &loc = ploc[pid];

        for (size_t i = 0; i < pf.size() && !found; i++) {
            E_Int fid = pf[i];
            const auto &pn = Fc[fid];
            const E_Float *fc = &fcenters[3*fid];

            if (pid == 276) {
                point_write("fc.im", fc[0], fc[1], fc[2]);
            }

            for (size_t j = 0; j < pn.size(); j++) {
                E_Int p = pn[j];
                E_Int q = pn[(j+1)%pn.size()];

                if (pid == 276) {
                    point_write("p.im", X[p], Y[p], Z[p]);
                    point_write("q.im", X[q], Y[q], Z[q]);
                }

                E_Float u, v, w;

                found = Triangle::is_point_inside(x, y, z,
                    X[p], Y[p], Z[p],
                    X[q], Y[q], Z[q],
                    fc[0], fc[1], fc[2],
                    u, v, w
                );

                if (found) {
                    loc.fid = fid;
                    loc.sub = j;
                    loc.bcrd[0] = u;
                    loc.bcrd[1] = v;
                    loc.bcrd[2] = w;

                    /*
                    if (pid == 107) {
                        printf("fid: %d - size: %d\n", fid, pn.size());
                        printf("u: %f - v: %f - w: %f\n", u, v, w);
                        printf("%f %f %f\n", fc[0], fc[1], fc[2]);
                    }
                    */

                    // on p
                    if      (Sign(1-u, NEAR_VERTEX_TOL) == 0) {
                        loc.v_idx = j;
                        loc.x = X[p]; loc.y = Y[p]; loc.z = Z[p];
                        loc.bcrd[0] = 1; loc.bcrd[1] = 0; loc.bcrd[2] = 0;
                    }
                    // on q
                    else if (Sign(1-v, NEAR_VERTEX_TOL) == 0) {
                        loc.v_idx = (j+1)%pn.size();
                        loc.x = X[q]; loc.y = Y[q]; loc.z = Z[q];
                        loc.bcrd[0] = 0; loc.bcrd[1] = 1; loc.bcrd[2] = 0;
                    }
                    // on edge {p, q}
                    else if (Sign(w, NEAR_EDGE_TOL) == 0) {
                        loc.e_idx = j;
                        loc.x = u*X[p] + (1-u)*X[q];
                        loc.y = u*Y[p] + (1-u)*Y[q];
                        loc.z = u*Z[p] + (1-u)*Z[q];
                        loc.bcrd[1] = 1-u; loc.bcrd[2] = 0;
                    }

                    break;
                }
            }
        }

        if (!found) {
            point_write("lost", x, y, z);
            write_ngon("bin", pf);
        }

        assert(found);
    }

    return ploc;
}

/*
void Smesh::correct_near_points_and_edges(Smesh &Sf,
    std::vector<PointLoc> &plocs)
{
    E_Int on_vertex = 0, on_edge = 0;
    for (size_t i = 0; i < plocs.size(); i++) {
        auto &ploc = plocs[i];

        E_Int fid = ploc.fid;
        assert(fid != -1);
        const auto &pn = Fc[fid];

        if (ploc.v_idx != -1) {
            on_vertex++;
            E_Int p = pn[ploc.v_idx];
            E_Float dx = X[p]-Sf.X[i];
            E_Float dy = Y[p]-Sf.Y[i];
            E_Float dz = Z[p]-Sf.Z[i];
            E_Float dist = dx*dx + dy*dy + dz*dz;
            if (dist >= Sf.min_pdist_squared) {
                fprintf(stderr, "Tight near-vertex situation!\n");
                point_write("mpoint", X[p], Y[p], Z[p]);
                point_write("spoint", Sf.X[i], Sf.Y[i], Sf.Z[i]);
                assert(0);
            } else {
                Sf.X[i] = X[p];
                Sf.Y[i] = Y[p];
                Sf.Z[i] = Z[p];
            }
        } else if (ploc.e_idx != -1) {
            on_edge++;
            E_Float u = ploc.bcrd[0];
            E_Float v = ploc.bcrd[1];
            E_Float w = ploc.bcrd[2];
            assert(Sign(w, NEAR_EDGE_TOL) == 0);
            v = 1 - u, w = 0;

            E_Int p = pn[ploc.e_idx];
            E_Int q = pn[(ploc.e_idx+1)%pn.size()];
            E_Float new_x = u*X[p] + v*X[q];
            E_Float new_y = u*Y[p] + v*Y[q];
            E_Float new_z = u*Z[p] + v*Z[q];

            E_Float dx = new_x-Sf.X[i];
            E_Float dy = new_y-Sf.Y[i];
            E_Float dz = new_z-Sf.Z[i];
            E_Float dist = dx*dx + dy*dy + dz*dz;
            if (dist >= Sf.min_pdist_squared) {
                fprintf(stderr, "Tight near-edge situation!\n");
                point_write("mpoint", X[p], Y[p], Z[p]);
                point_write("spoint", Sf.X[i], Sf.Y[i], Sf.Z[i]);
                assert(0);
            } else {
                ploc.bcrd[1] = v; ploc.bcrd[2] = 0;

                Sf.X[i] = new_x;
                Sf.Y[i] = new_y;
                Sf.Z[i] = new_z;
            }
        }
    }
    printf("on vertex: %d - on edge: %d\n", on_vertex, on_edge);
}
*/

void Smesh::make_bbox()
{
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

    E_Float dx = xmax - xmin;
    E_Float dy = ymax - ymin;
    E_Float dz = zmax - zmin;

    //xmin = xmin - dx*0.01;
    //ymin = ymin - dy*0.01;
    //zmin = zmin - dz*0.01;
    //xmax = xmax + dx*0.01;
    //ymax = ymax + dy*0.01;
    //zmax = zmax + dz*0.01;

    //HX = (dx != 0) ? (xmax - xmin) / NX : 1;
    //HY = (dy != 0) ? (ymax - ymin) / NY : 1;
    //HZ = (dz != 0) ? (zmax - zmin) / NZ : 1;

    box.xmin = xmin, box.ymin = ymin, box.zmin = zmin;
    box.xmax = xmax, box.ymax = ymax, box.zmax = zmax;
}

inline
void Smesh::bin_face(E_Int fid)
{
    const auto &pn = Fc[fid];

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

void Smesh::hash_faces()
{
    bin_faces.clear();

    for (E_Int fid = 0; fid < nf; fid++) {
        bin_face(fid);
    }
}
