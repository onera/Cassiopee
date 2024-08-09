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

struct IMesh;

struct o_edge {
    Int p, q;

    o_edge(Int P, Int Q);
};

struct Smesh {
    Int np, ne, nf;

    std::vector<Float> X, Y, Z;
    std::vector<std::vector<Int>> P2F;
    std::vector<std::vector<Int>> P2E;
    
    std::vector<o_edge> E;
    std::vector<std::array<Int, 2>> E2F;

    std::vector<std::vector<Int>> F;
    std::vector<std::vector<Int>> F2E;
    std::vector<std::vector<Int>> F2F;

    Int M_np;
    Int M_ne;
    Int M_nf;

    std::map<Int, Int> g2lp;
    std::map<Int, Int> l2gp;

    std::map<Int, Int> g2lf;
    std::map<Int, Int> l2gf;

    std::map<Int, Int> g2le;
    std::map<Int, Int> l2ge;

    // Adaptation
    void get_leaves(Int face, std::vector<Int> &leaves) const;


    std::map<Int, std::vector<Int>> fchildren;
    std::set<Int> factive;
    std::vector<Int> flevel;

    std::map<Int, std::vector<Int>> echildren;
    std::set<Int> eactive;
    std::vector<Int> elevel;

    Smesh() {};

    Smesh(const char *fname);

    Smesh(const IMesh &M);
    
    Smesh(const IMesh &M, const std::vector<Int> &faces);

    bool ccw_oriented(Int face);

    void make_edges();

    void make_point_faces();
    
    void make_point_edges();

    inline bool edge_is_active(Int edge) const
    { return eactive.find(edge) != eactive.end(); }

    inline bool face_is_active(Int face) const
    { return factive.find(face) != factive.end(); }

    size_t refine(Smesh &M);

    std::vector<pointFace> locate(Float x, Float y, Float z) const;

    void write_faces(const char *fname, const std::vector<Int> &faces);

    void write_ngon(const char *fname);

    inline bool face_is_tri(Int fid) const { return F[fid].size() == 3; }

    inline bool face_is_quad(Int fid) const { return F[fid].size() == 4; }

    bool face_contains_Mface(Int face, Int mface, const Smesh &M) const;

    Int face_contains_point(Int face, Float x, Float y, Float z) const;

    std::vector<Int> smooth_ref_data(std::map<Int, std::vector<Int>> &sensor);

    std::vector<Int> prepare_for_refinement(const std::vector<Int> &ref_data);

    void refine_faces(const std::vector<Int> &ref_faces);

    std::vector<Int> get_active_neighbours(Int face);

    inline Int get_neighbour(Int face, Int edge) const
    {
        assert(E2F[edge][0] == face || E2F[edge][1] == face);
        return (E2F[edge][0] == face) ? E2F[edge][1] : E2F[edge][0];
    }

    void resize_point_data(size_t nref_faces);

    void resize_edge_data(size_t nref_faces);

    void resize_face_data(size_t nref_faces);

    void refine_tri(Int tri);

    void refine_quad(Int quad);

    void refine_edge(Int edge);
    
    inline Int get_edge_center(Int edge)
    {
        assert(echildren.find(edge) != echildren.end());
        return E[echildren[edge][0]].q;
    }

    void init_adaptation_data();

    Smesh extract_conformized();

    bool faces_are_dups(Int face, Int mface, const Smesh &M);
};
