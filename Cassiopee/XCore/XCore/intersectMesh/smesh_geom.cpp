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
#include "io.h"
#include "ray.h"

struct HitData {
    std::set<E_Int> tested;
    std::set<E_Int> fids;
    std::set<Point3D> hits;
};

void Smesh::intersect_ray(const Ray &ray, E_Int node_idx, HitData &hit_data) const
{
    const BVH_node &node = bvh_nodes[node_idx];
    if (!ray.intersect_AABB(node.box)) return;
    if (!node.is_leaf()) {
        intersect_ray(ray, node.left_node, hit_data);
        intersect_ray(ray, node.left_node+1, hit_data);
        return;
    }
    for (E_Int i = 0; i < node.tri_count; i++) {
        E_Int tri = tri_idx[node.first_tri_idx+i];
        hit_data.tested.insert(tri);
        const auto &pn = F[tri];
        //write_face("tri.im", tri);
        E_Float u, v, w, t, x, y, z;
        u = v = w = 2.0;
        bool hit = ray.intersect_triangle({X[pn[0]], Y[pn[0]], Z[pn[0]]},
            {X[pn[1]], Y[pn[1]], Z[pn[1]]}, {X[pn[2]], Y[pn[2]], Z[pn[2]]},
            u, v, w, t, x, y, z);
        
        if (hit) {
            hit_data.fids.insert(tri);
            hit_data.hits.insert({x, y, z});
        }
    }
}

bool Smesh::is_point_inside(E_Float px, E_Float py, E_Float pz) const
{
    Ray ray(px, py, pz, 1, 0, 0);
    return is_point_inside(ray);
}

bool Smesh::is_point_inside(const Ray &ray) const
{
    //if (!box.is_point_inside(ray.o)) return false;
    HitData hit_data;
    intersect_ray(ray, 0, hit_data);
    return hit_data.hits.size() % 2 == 1;
}

void Smesh::replace_by_projections(const std::vector<E_Int> &pids,
    const std::vector<PointLoc> &plocs)
{
    for (size_t i = 0; i < pids.size(); i++) {
        E_Int pid = pids[i];
        const auto &ploc = plocs[pid];
        if (ploc.v_idx != -1 || ploc.e_idx != -1) {
            E_Float dx = X[pid]-ploc.x;
            E_Float dy = Y[pid]-ploc.y;
            E_Float dz = Z[pid]-ploc.z;
            E_Float d = sqrt(dx*dx + dy*dy + dz*dz);
            if (d >= min_pdist/10) {
                fprintf(stderr, "Tight near vertex/edge situation!\n");
                point_write("orig.im", X[pid], Y[pid], Z[pid]);
                point_write("dest.im", ploc.x, ploc.y, ploc.z);
                assert(0);
            }
            X[pid] = ploc.x;
            Y[pid] = ploc.y;
            Z[pid] = ploc.z;
        }
    }
}

void Smesh::compute_min_distance_between_points()
{
    min_pdist = EFLOATMAX;
    size_t ndists = 0;

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

            E_Float dist = sqrt(dx*dx + dy*dy + dz*dz);

            if (dist < min_pdist) min_pdist = dist;

            ndists++;
        }
    }

    assert(ndists == (size_t)np*((size_t)np-1)/2);
}

void Smesh::get_unit_projected_direction(E_Int fid, const E_Float D[3],
    E_Float proj[3]) const
{
    assert(fid >= 0);
    assert(fid < nf);

    // Unit normal
    const E_Float *fN = &fnormals[3*fid];
    assert(Sign(K_MATH::norm(fN, 3)) != 0);

    E_Float dp = K_MATH::dot(D, fN, 3); 

    proj[0] = D[0] - dp * fN[0];
    proj[1] = D[1] - dp * fN[1];
    proj[2] = D[2] - dp * fN[2];
    E_Float NORM = K_MATH::norm(proj, 3);
    proj[0] /= NORM, proj[1] /= NORM, proj[2] /= NORM;
}

E_Int Smesh::deduce_face(const std::vector<E_Int> &pf,
    E_Float ox, E_Float oy, E_Float oz, E_Float D[3], 
    E_Int last_vertex, E_Int last_edge, E_Int eid) const
{
    // Intersect the projection of D with all the faces in pf
    // At least one intersection must exist
    // Return the face with the earliest intersection

    // For debugging
    E_Int faces_hit = 0;

    E_Float t_min = EFLOATMAX;
    E_Int ret_face = -1;

    for (auto fid : pf) {

        // Compute the unit projection of D on this face

        E_Float proj[3];
        get_unit_projected_direction(fid, D, proj);

        const auto &pn = Fc[fid];
        const auto &pe = F2E[fid];
        assert(pn.size() == pe.size());

        for (size_t i = 0; i < pn.size(); i++) {

            E_Int p = pn[i];
            E_Int q = pn[(i+1)%pn.size()];
            E_Int e = pe[i];

            if (p == last_vertex || q == last_vertex || e == last_edge)
                continue;

            E_Float t, s;

            //if (eid == 133) {
            //    point_write("p.im", X[p], Y[p], Z[p]);
            //    point_write("q.im", X[q], Y[q], Z[q]);
            //}

            bool hit = ray_edge_intersect(ox, oy, oz,
                proj[0], proj[1], proj[2],
                X[p], Y[p], Z[p],
                X[q], Y[q], Z[q],
                t, s
            );

            if (hit) {
                faces_hit += 1;

                if (t < t_min) {
                    t_min = t;
                    ret_face = fid;
                }

                // Hit an edge of the face, stop
                break;
            }
        }
    }

    // We must have hit a face
    if (faces_hit == 0) {
        write_ngon("pf.im", pf);
        edge_write("ray.im", ox, oy, oz, ox+D[0], oy+D[1], oz+D[2]);
    }
    assert(faces_hit > 0);

    return ret_face;
}
