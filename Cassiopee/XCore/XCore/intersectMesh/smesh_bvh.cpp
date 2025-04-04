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
#include "BVH.h"
#include "ray.h"
#include "primitives.h"

void Smesh::make_BVH()
{
    bvh_nodes.clear();
    bvh_nodes.resize(2*nf - 1);
    root_node_idx = 0, nodes_used = 1;

    tri_idx.clear();
    tri_idx.resize(nf);
    for (E_Int i = 0; i < nf; i++) tri_idx[i] = i;
    assert(fcenters.size() == (size_t)nf*3);

    // Assign all triangles to root node
    BVH_node &root = bvh_nodes[root_node_idx];
    root.left_node = 1;
    root.first_tri_idx = 0;
    root.tri_count = nf;
    update_node_bounds(root_node_idx);
    BVH_subdivide(root_node_idx);
}

void Smesh::update_node_bounds(E_Int node_idx)
{
    BVH_node &node = bvh_nodes[node_idx];
    node.box.xmin = node.box.ymin = node.box.zmin = EFLOATMAX;
    node.box.xmax = node.box.ymax = node.box.zmax = EFLOATMIN;
    for (E_Int first = node.first_tri_idx, i = 0; i < node.tri_count; i++) {
        E_Int leaf_tri_idx = tri_idx[first + i];
        const auto &pn = F[leaf_tri_idx];
        for (E_Int p : pn) {
            node.box.xmin = std::min(node.box.xmin, X[p]);
            node.box.ymin = std::min(node.box.ymin, Y[p]);
            node.box.zmin = std::min(node.box.zmin, Z[p]);
            node.box.xmax = std::max(node.box.xmax, X[p]);
            node.box.ymax = std::max(node.box.ymax, Y[p]);
            node.box.zmax = std::max(node.box.zmax, Z[p]);
        }
        node.box.dx = node.box.xmax - node.box.xmin;
        node.box.dy = node.box.ymax - node.box.ymin;
        node.box.dz = node.box.zmax - node.box.zmin;
    }
}

void Smesh::BVH_subdivide(E_Int node_idx)
{
    // Terminate recursion
    BVH_node &node = bvh_nodes[node_idx];
    if (node.tri_count <= MAX_TRIS_PER_BVH_LEAF) return;

    // Determine split axis and position
    int axis;
    E_Float split_pos;
    if (node.box.dx > node.box.dy && node.box.dx > node.box.dz) {
        axis = 0;
        split_pos = node.box.xmin + node.box.dx*0.5;
    }
    else if (node.box.dy > node.box.dz) {
        axis = 1;
        split_pos = node.box.ymin + node.box.dy*0.5;
    }
    else {
        axis = 2;
        split_pos = node.box.zmin + node.box.dz*0.5;
    }

    // In-place partition
    E_Int i = node.first_tri_idx;
    E_Int j = i + node.tri_count - 1;
    while (i <= j) {
        const E_Float *fc = &fcenters[3*tri_idx[i]];
        if (fc[axis] < split_pos)
            i++;
        else
            std::swap(tri_idx[i], tri_idx[j--]);
    }

    // Abort split if one of the sides is empty
    E_Int left_count = i - node.first_tri_idx;
    if (left_count == 0 || left_count == node.tri_count) return;

    // Create child nodes
    E_Int left_child_idx = nodes_used++;
    E_Int right_child_idx = nodes_used++;
    bvh_nodes[left_child_idx].first_tri_idx = node.first_tri_idx;
    bvh_nodes[left_child_idx].tri_count = left_count;
    bvh_nodes[right_child_idx].first_tri_idx = i;
    bvh_nodes[right_child_idx].tri_count = node.tri_count - left_count;
    node.left_node = left_child_idx;
    node.tri_count = 0;
    update_node_bounds(left_child_idx);
    update_node_bounds(right_child_idx);

    // Recurse
    BVH_subdivide(left_child_idx);
    BVH_subdivide(right_child_idx);
}

void Smesh::make_BVH(const std::set<E_Int> &fids)
{
    bvh_nodes.clear();
    size_t NF = fids.size();
    bvh_nodes.resize(2*NF - 1);
    root_node_idx = 0, nodes_used = 1;

    tri_idx.clear();
    tri_idx.reserve(NF);
    for (E_Int fid : fids) tri_idx.push_back(fid);

    // Assign all triangles to root node
    BVH_node &root = bvh_nodes[root_node_idx];
    root.left_node = 1;
    root.first_tri_idx = 0;
    root.tri_count = NF;
    update_node_bounds(root_node_idx);
    BVH_subdivide(root_node_idx);
}

static
bool ray_intersect_AABB(E_Float ox, E_Float oy, E_Float oz,
    E_Float dx, E_Float dy, E_Float dz, const AABB &box)
{
    E_Float tmin = 0;
    E_Float tmax = EFLOATMAX;

    E_Float origin[3] = { ox, oy, oz };
    E_Float direction[3] = { dx, dy, dz };
    E_Float boxMin[3] = { box.xmin, box.ymin, box.zmin };
    E_Float boxMax[3] = { box.xmax, box.ymax, box.zmax };

    for (int i = 0; i < 3; i++) {
        E_Float o = origin[i];
        E_Float d = direction[i];
        E_Float bmin = boxMin[i];
        E_Float bmax = boxMax[i];

        if (d != 0) {
            E_Float t1 = (bmin - o) / d;
            E_Float t2 = (bmax - o) / d;

            if (t1 > t2) { E_Float temp = t1; t1 = t2; t2 = temp; }

            tmin = (t1 > tmin) ? t1 : tmin;
            tmax = (t2 < tmax) ? t2 : tmax;

            if (tmin > tmax) return false;  // No intersection
        } else {
            if (o < bmin || o > bmax) return false;  // Parallel and outside slab
        }
    }

    return true;
}

/*
void Smesh::ray_intersect_BVH(E_Float ox, E_Float oy, E_Float oz,
    E_Float dx, E_Float dy, E_Float dz, E_Int node_idx,
    std::vector<PointLoc> &plocs) const
{
    const BVH_node &node = bvh_nodes[node_idx];
    if (!ray_intersect_AABB(ox, oy, oz, dx, dy, dz, node.box)) return;
    if (node.is_leaf()) {
        for (E_Int i = 0; i < node.tri_count; i++) {
            E_Int tri = tri_idx[node.first_tri_idx+i];
            const auto &pn = Fc[tri];
            const E_Float *fc = &fcenters[3*tri];

            for (size_t j = 0; j < pn.size(); j++) {
                E_Int p = pn[j];
                E_Int q = pn[(j+1)%pn.size()];
                E_Float u, v, w, t, x, y, z;

                bool hit = MollerTrumboreAnyDir(
                    ox, oy, oz, dx, dy, dz,
                    X[p], Y[p], Z[p],
                    X[q], Y[q], Z[q],
                    fc[0], fc[1], fc[2],
                    u, v, w, t, x, y, z
                );
                
                if (hit) {
                    PointLoc ploc;
                    ploc.fid = tri;
                    ploc.sub = j;
                    ploc.bcrd[0] = u;
                    ploc.bcrd[1] = v;
                    ploc.bcrd[2] = w;
                    ploc.t = t;
                    ploc.x = x;
                    ploc.y = y;
                    ploc.z = z;

                    // on p
                    if      (Sign(1-u, NEAR_VERTEX_TOL) == 0) {
                        ploc.v_idx = j;
                        ploc.bcrd[0] = 1, ploc.bcrd[1] = 0, ploc.bcrd[2] = 0;
                        ploc.x = X[p];
                        ploc.y = Y[p];
                        ploc.z = Z[p];
                    }
                    // on q
                    else if (Sign(1-v, NEAR_VERTEX_TOL) == 0) {
                        ploc.v_idx = (j+1)%pn.size();
                        ploc.bcrd[0] = 0, ploc.bcrd[1] = 1, ploc.bcrd[2] = 0;
                        ploc.x = X[q];
                        ploc.y = Y[q];
                        ploc.z = Z[q];
                    }
                    // on edge {p, q}
                    else if (Sign(w, NEAR_EDGE_TOL) == 0) {
                        ploc.e_idx = j;
                        ploc.bcrd[0] = u, ploc.bcrd[1] = 1-u, ploc.bcrd[2] = 0;
                        ploc.x = u*X[p] + (1-u)*X[q];
                        ploc.y = u*Y[p] + (1-u)*Y[q];
                        ploc.z = u*Z[p] + (1-u)*Z[q];
                    }

                    plocs.push_back(ploc);
                    
                    break;
                }
            }
        }
        return ;
    }

    ray_intersect_BVH(ox, oy, oz, dx, dy, dz, node.left_node, plocs);
    ray_intersect_BVH(ox, oy, oz, dx, dy, dz, node.left_node+1, plocs);
}
*/

