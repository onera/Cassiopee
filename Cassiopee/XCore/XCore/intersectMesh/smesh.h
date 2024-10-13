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
#pragma once

#include <vector>
#include <array>
#include <map>
#include <set>
#include <cassert>

#include "point.h"
#include "xcore.h"
#include "AABB.h"

struct IMesh;

struct o_edge {
    E_Int p, q;

    o_edge(E_Int P, E_Int Q);
};

struct u_edge {
    E_Int p, q;

    u_edge(E_Int P, E_Int Q)
    : p(std::min(P, Q)), q(std::max(P, Q))
    {}

    bool operator<(const u_edge& e) const
    {
        return (p < e.p) || (p == e.p && q < e.q);
    }
};

struct Smesh {
    E_Int np, ne, nf;

    std::vector<E_Float> X, Y, Z;
    std::vector<std::vector<E_Int>> P2F;
    std::vector<std::vector<E_Int>> P2E;
    
    std::vector<o_edge> E;
    std::vector<std::array<E_Int, 2>> E2F;

    std::vector<std::vector<E_Int>> F;
    std::vector<std::vector<E_Int>> F2E;
    std::vector<std::vector<E_Int>> F2F;

    std::map<E_Int, E_Int> g2lp;
    std::map<E_Int, E_Int> l2gp;

    std::map<E_Int, E_Int> g2lf;
    std::map<E_Int, E_Int> l2gf;

    std::vector<E_Float> fnormals;
    std::vector<E_Float> pnormals;

    E_Int NX, NY, NZ, NXY, NXYZ;
    E_Float xmin, xmax, ymin, ymax, zmin, zmax;
    E_Float HX, HY, HZ;

    std::vector<std::vector<E_Int>> bin_faces;

    E_Float min_pdist_squared = EFLOATMIN;

    E_Float NEAR_VERTEX_TOL = 1e-3;
    E_Float NEAR_EDGE_TOL = 1e-3;

    std::map<u_edge, E_Int> ecenter;

    Smesh();

    void resize_for_refinement(size_t nref_faces);

    void refine(const std::map<E_Int, std::vector<PointData>> &sensor);

    std::set<E_Int> extract_bounding_faces(const Smesh &Sf,
        const std::vector<PointLoc> &plocs) const;

    void correct_near_points_and_edges(Smesh &Sf, std::vector<PointLoc> &plocs);
    
    E_Int deduce_face(const std::vector<E_Int> &pf,
        E_Float ox, E_Float oy, E_Float oz, E_Float D[3], 
        E_Int last_vertex, E_Int last_edge) const;

    void get_unit_projected_direction(E_Int fid, const E_Float D[3],
        E_Float proj[3]) const;
    
    void get_shared_faces(const PointLoc &loc, std::vector<E_Int> &ret,
        E_Int &pid, E_Int &eid) const;
    
    void make_fnormals();
    void make_pnormals();

    void make_bbox();
    void hash_faces();

    void compute_min_distance_between_points();

    void get_leaves(E_Int face, std::vector<E_Int> &leaves) const;

    std::map<E_Int, std::vector<std::array<E_Int, 3>>> fchildren;

    Smesh extract_smesh(const std::set<E_Int> &fids, bool is_planar=true);

    Smesh(const IMesh &M, bool is_planar=true);
    
    Smesh(const IMesh &M, const std::vector<E_Int> &faces);

    void make_edges(bool is_planar);

    void make_point_faces();
    
    void make_point_edges();

    std::vector<PointLoc> locate(const Smesh &Sf) const;

    void write_points(const char *fname, const std::set<E_Int> &pset) const;
    
    void write_points(const char *fname, const std::vector<E_Int> &pids) const;

    void write_edges(const char *fname, const std::set<E_Int> &eids) const;

    void write_ngon(const char *fname, const std::set<E_Int> &fset) const;

    void write_ngon(const char *fname, const std::vector<E_Int> &faces) const;

    void write_ngon(const char *fname);

    void refine_faces(const std::vector<E_Int> &ref_faces);

    void refine_tri(E_Int fid);

    void refine_edge(const u_edge &e);
    
    Smesh extract_conformized();

    inline E_Int get_voxel(E_Int i, E_Int j, E_Int k) const
    {
        return i + NX*j + NXY*k;
    }
};

