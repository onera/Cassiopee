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

    E_Int M_np;
    E_Int M_ne;
    E_Int M_nf;

    std::map<E_Int, E_Int> g2lp;
    std::map<E_Int, E_Int> l2gp;

    std::map<E_Int, E_Int> g2lf;
    std::map<E_Int, E_Int> l2gf;

    std::map<E_Int, E_Int> g2le;
    std::map<E_Int, E_Int> l2ge;

    std::vector<E_Float> fnormals;
    std::vector<E_Float> pnormals;

    E_Int NX, NY, NZ, NXY;
    E_Float xmin, xmax, ymin, ymax, zmin, zmax;
    E_Float HX, HY, HZ;

    std::map<E_Int, std::vector<E_Int>> fmap;

    void make_fnormals();
    void make_pnormals();

    void make_bbox();
    void hash_faces();

    AABB AABB_face(const std::vector<E_Int> &pn) const;

    // Adaptation
    void get_leaves(E_Int face, std::vector<E_Int> &leaves) const;


    std::map<E_Int, std::vector<E_Int>> fchildren;
    std::set<E_Int> factive;
    std::vector<E_Int> flevel;

    std::map<E_Int, std::vector<E_Int>> echildren;
    std::set<E_Int> eactive;
    std::vector<E_Int> elevel;

    Smesh() {};

    Smesh(const char *fname);

    Smesh(const IMesh &M);
    
    Smesh(const IMesh &M, const std::vector<E_Int> &faces);

    bool ccw_oriented(E_Int face);

    void make_edges();

    void make_point_faces();
    
    void make_point_faces_all();
    
    void make_point_edges();

    inline bool edge_is_active(E_Int edge) const
    { return eactive.find(edge) != eactive.end(); }

    inline bool face_is_active(E_Int face) const
    { return factive.find(face) != factive.end(); }

    //size_t refine(Smesh &M);

    //std::vector<pointFace> locate(E_Float x, E_Float y, E_Float z) const;

    void write_faces(const char *fname, const std::vector<E_Int> &faces) const;

    void write_ngon(const char *fname);

    inline bool face_is_tri(E_Int fid) const { return F[fid].size() == 3; }

    inline bool face_is_quad(E_Int fid) const { return F[fid].size() == 4; }

    bool face_contains_Mface(E_Int face, E_Int mface, const Smesh &M) const;

    E_Int face_contains_point(E_Int face, E_Float x, E_Float y, E_Float z) const;

    std::vector<E_Int> smooth_ref_data(std::map<E_Int, std::vector<E_Int>> &sensor);

    std::vector<E_Int> prepare_for_refinement(const std::vector<E_Int> &ref_data);

    void refine_faces(const std::vector<E_Int> &ref_faces);

    std::vector<E_Int> get_active_neighbours(E_Int face);

    inline E_Int get_neighbour(E_Int face, E_Int edge) const
    {
        assert(E2F[edge][0] == face || E2F[edge][1] == face);
        return (E2F[edge][0] == face) ? E2F[edge][1] : E2F[edge][0];
    }

    void resize_point_data(size_t nref_faces);

    void resize_edge_data(size_t nref_faces);

    void resize_face_data(size_t nref_faces);

    void refine_tri(E_Int tri);

    void refine_quad(E_Int quad);

    void refine_edge(E_Int edge);
    
    inline E_Int get_edge_center(E_Int edge)
    {
        assert(echildren.find(edge) != echildren.end());
        return E[echildren[edge][0]].q;
    }

    void init_adaptation_data();

    Smesh extract_conformized();

    bool faces_are_dups(E_Int face, E_Int mface, const Smesh &M);
};
