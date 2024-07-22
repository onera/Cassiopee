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

#include "point.h"
#include "xcore.h"
#include "common/common.h"

#define OUT 0
#define IN 1

struct Ray;
struct Smesh;

struct UEdge {
    Int p, q;

    UEdge(Int P, Int Q);

    bool operator<(const UEdge &E) const;
};

struct IMesh {
    Int np, ne, nf, nc;

    std::vector<Float> X, Y, Z;
    std::vector<std::vector<Int>> P2F;

    std::vector<std::array<Int, 2>> E;
    
    std::vector<std::vector<Int>> F;
    std::vector<std::vector<Int>> F2E;

    std::vector<std::vector<Int>> C;

    std::vector<Int> skin;

    std::set<Int> patch;

    Float xmin, ymin, zmin;
    Float xmax, ymax, zmax;
    Float dmin, dmax;
    Float DX, DY, DZ;
    Int NBIN;

    std::map<Int, std::set<Int>> bin_faces;

    IMesh();

    IMesh(const char *fname);

    IMesh(K_FLD::FldArrayI &cn, Float *X, Float *Y, Float *Z, Int npts);

    void make_skin();

    void make_bbox();

    void hash_skin();

    void make_point_faces();

    void make_edges();

    inline bool face_is_quad(Int face) const { return F[face].size() == 4; }
    
    inline bool face_is_tri(Int face) const { return F[face].size() == 3; }

    void write_ngon(const char *fname);

    void write_faces(const char *fname, const std::vector<Int> &faces);

    bool is_point_inside(Float px, Float py, Float pz);

    IMesh reconstruct_after_smesh_adaptation(const Smesh &Mf, Int patchc);

    // Adaptation
    void init_adaptation_data();

    std::vector<Int> smooth_ref_data(
        const std::map<Int, std::vector<Int>> &sensor);

    std::vector<Int> prepare_for_refinement(
        const std::vector<Int> &ref_data);

    std::set<Int> factive;
    std::map<Int, std::vector<Int>> fchildren;
    std::vector<Int> flevel;

    std::map<UEdge, Int> ecenter;

    size_t refine(std::set<Int> &mpatch, std::set<Int> &spatch, IMesh &S);

    std::vector<pointFace> locate(Float x, Float y,
        const std::set<Int> &patch) const;
    
    inline bool face_is_active(Int face) const
    { return factive.find(face) != factive.end(); }

    bool faces_are_dups(Int mface, Int sface, const IMesh &S);

    void refine_faces(const std::vector<Int> &ref_faces);

    bool face_contains_sface(Int face, Int sface, const IMesh &S) const;

    void resize_point_data(size_t nref_faces);

    void resize_face_data(size_t nref_faces);

    void refine_quad(Int quad);

    void refine_tri(Int tri);

    void refine_edge(const UEdge &edge);

    Int face_contains_point(Int face, Float x, Float y) const;

    IMesh extract_conformized();

    void get_fleaves(Int face, std::vector<Int> &fleaves);

    PyObject *export_karray();

    /* TOPO */

    Int orient_skin(Int normal_direction);

    void flag_and_get_external_faces(std::vector<Int> &fflags,
        std::vector<Int> &efaces);
    
    void extract_nface_of_kept_pgs(const std::vector<bool> &kept_pgs,
        std::vector<Int> &NFACE, std::vector<Int> &cxadj,
        std::vector<Int> &cells);
    
    void flag_marked_external_cells(const std::vector<Int> &cells,
        const std::vector<Int> &fflags, std::vector<Int> &cflags);
    
    void flag_all_external_cells(const std::vector<Int> &fflags,
        std::vector<Int> &cflags);
    
    Int orient_boundary(Int ncells, Int *efadj, Int *exadj, Int nefaces,
        Int *fneis, Int *efaces, std::vector<Int> &forient,
        const std::vector<Int> &cflags, const std::vector<Int> &fflags,
        Int *cells, Int normal_direction);
    
    void compute_cell_volume(Int cell, Float &vol, Int refIdx);
};

Int meshes_mutual_refinement(IMesh &M, IMesh &S);