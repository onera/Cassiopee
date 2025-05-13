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
#pragma once

#include <vector>
#include <array>
#include <map>
#include <set>

#include "xcore.h"
#include "u_edge.h"
#include "point.h"
#include "BVH.h"
#include "AABB.h"

#define LEFT 0
#define RIGHT 1

struct IMesh;
struct Ray;
struct HitData;

struct o_edge {
    E_Int p, q;

    o_edge(E_Int P, E_Int Q);
};

struct Smesh {

    // Universal data

    E_Int np, nf;
    std::vector<E_Float> X, Y, Z;
    std::vector<std::vector<E_Int>> F;

    // Conformal data
    E_Int ne;
    std::vector<std::vector<E_Int>> Fc;
    std::vector<std::vector<E_Int>> P2F;
    std::vector<std::vector<E_Int>> P2E;
    std::vector<o_edge> E;
    std::vector<std::array<E_Int, 2>> E2F;
    std::vector<std::vector<E_Int>> F2E;
    std::vector<std::vector<E_Int>> F2F;

    void clear_conformal_data();

    // Link to parent mesh

    std::map<E_Int, E_Int> g2lp;
    std::map<E_Int, E_Int> l2gp;
    std::map<E_Int, E_Int> g2lf;
    std::map<E_Int, E_Int> l2gf;

    void reconstruct(IMesh &M);

    // Constructors

    Smesh();
    Smesh(const Smesh &Mf, const std::set<E_Int> &fids, bool check_Euler=true);
    Smesh(const IMesh &M, const std::vector<E_Int> &fids, bool check_Euler=true);
    static Smesh Smesh_from_point_tags(const IMesh &M, const E_Float *ptag,
        bool check_Euler=true);
    static Smesh Smesh_from_mesh_skin(const IMesh &M,
        const std::vector<E_Int> &skin, bool check_Euler=true);
    static Smesh Smesh_from_mesh_patch(const IMesh &M,
        bool check_Euler=true);
    static Smesh Smesh_from_tagged_faces(const IMesh &M, bool check_Euler);
    static Smesh make_sub(const Smesh &Mf, const std::set<E_Int> &bfaces,
        bool check_Euler);
    void make_edges();
    void clear();
    void tag_faces(IMesh &M) const;

    // Geometry

    std::vector<E_Float> fcenters;
    std::vector<E_Float> fnormals;
    std::vector<E_Float> pnormals;
    E_Float min_pdist = EFLOATMIN;
    E_Float NEAR_VERTEX_TOL = 1e-3;
    E_Float NEAR_EDGE_TOL = 1e-3;

    bool is_point_in_3D_polygon(E_Float x, E_Float y, E_Float z, E_Int fid) const;
    bool is_point_on_a_polygon_edge(E_Float x, E_Float y, E_Float z, E_Int fid,
        PointLoc &ploc, E_Float min_pdist) const;
    bool is_point_a_polygon_vertex(E_Float x, E_Float y, E_Float z, E_Int fid,
        PointLoc &ploc, E_Float min_pdist) const;
    bool is_point_inside(E_Float px, E_Float py, E_Float pz) const;
    bool is_point_inside(const Ray &ray) const;
    void intersect_ray(const Ray &ray, E_Int node_idx, HitData &hit_data) const;
    
    void make_fcenters();
    void make_fnormals();
    void make_pnormals();
    void make_point_faces();
    void make_point_edges();
    std::vector<PointLoc> locate(const Smesh &Sf) const;
    std::vector<PointLoc> locate2(const Smesh &Sf) const;
    void correct_near_points_and_edges(Smesh &Sf, std::vector<PointLoc> &plocs);
    void get_unit_projected_direction(E_Int fid, const E_Float D[3],
        E_Float proj[3]) const;
    void compute_min_distance_between_points();
    std::vector<PointLoc> project(const Smesh &Mf,
        const std::vector<E_Int> &mpids) const;
    void replace_by_projections(const std::vector<E_Int> &pids,
        const std::vector<PointLoc> &plocs);
    inline void compute_face_center(E_Int fid);

    // Hash
    AABB box;
    E_Int NX, NY, NZ, NXY, NXYZ;
    E_Float xmin, xmax, ymin, ymax, zmin, zmax;
    E_Float HX, HY, HZ;
    std::map<E_Int, std::vector<E_Int>> bin_faces;
    
    void make_bbox();
    inline void bin_face(E_Int fid);
    void hash_faces();
    inline E_Int get_voxel(E_Int i, E_Int j, E_Int k) const
    {
        return i + NX*j + NXY*k;
    }


    // BVH
        
    E_Int root_node_idx, nodes_used;
    std::vector<E_Int> tri_idx;
    std::vector<BVH_node> bvh_nodes;
    static const E_Int MAX_TRIS_PER_BVH_LEAF = 2;
    void make_BVH();
    void make_BVH(const std::set<E_Int> &fids);
    void update_node_bounds(E_Int node_idx);
    void BVH_subdivide(E_Int node_idx);
    void ray_intersect_BVH(E_Float ox, E_Float oy, E_Float oz,
        E_Float dx, E_Float dy, E_Float dz, E_Int node_idx,
        std::vector<PointLoc> &plocs) const;

    // Adaptation

    std::map<u_edge, E_Int> ecenter;
    std::map<E_Int, std::vector<std::vector<E_Int>>> fchildren;
    E_Int np_before_adapt;
    E_Int nf_before_adapt;
    
    void resize_for_refinement(size_t nref_faces);
    void refine(const std::vector<E_Int> &ref_faces);
    void get_leaves(E_Int face, std::vector<E_Int> &leaves) const;
    void refine_tri(E_Int fid);
    void refine_edge(const u_edge &e);
    void update_plocs(const std::vector<E_Int> &parents,
        std::vector<PointLoc> &plocs);
    void conformize();
    void get_edge_centers(E_Int p, E_Int q, std::vector<E_Int> &edge_centers);
    std::vector<E_Int> deduce_ref_faces(const std::vector<E_Int> &mpids,
        const std::vector<PointLoc> &plocs_m, const Smesh &Mf,
        std::vector<E_Int> &ref_faces);

    // Topology

    bool check_Euler = true;

    std::set<E_Int> extract_covering_faces(const Smesh &Sf,
        const std::vector<PointLoc> &plocs) const;
    E_Int deduce_face(const std::vector<E_Int> &pf,
        E_Float ox, E_Float oy, E_Float oz, E_Float D[3], 
        E_Int last_vertex, E_Int last_edge, E_Int eid) const;
    void get_shared_faces(const PointLoc &loc, std::vector<E_Int> &ret,
        E_Int &pid, E_Int &eid) const;
    Smesh extract_smesh(const std::set<E_Int> &fids, bool check_Euler=true);
    Smesh extract_conformized();

    // IO

    void write_points(const char *fname, const std::set<E_Int> &pset) const;
    void write_points(const char *fname, const std::vector<E_Int> &pids) const;
    void write_edges(const char *fname, const std::set<E_Int> &eids) const;
    void write_ngon(const char *fname, const std::set<E_Int> &fset) const;
    void write_ngon(const char *fname, const std::vector<E_Int> &faces) const;
    void write_ngon(const char *fname) const;
    void write_face(const char *fname, E_Int fid) const;
};

