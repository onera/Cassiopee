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

#define OUT 0
#define IN 1

struct Ray;
struct Smesh;

struct UEdge {
    E_Int p, q;

    UEdge(E_Int P, E_Int Q);

    bool operator<(const UEdge &E) const;
};

struct Mesh {
    E_Int np, ne, nf, nc;

    std::vector<E_Float> X, Y, Z;
    std::vector<std::vector<E_Int>> P2F;

    std::vector<std::array<E_Int, 2>> E;
    
    std::vector<std::vector<E_Int>> F;
    std::vector<std::vector<E_Int>> F2E;

    std::vector<std::vector<E_Int>> C;

    std::vector<E_Int> skin;

    std::set<E_Int> patch;

    E_Float xmin, ymin, zmin;
    E_Float xmax, ymax, zmax;
    E_Float dmin, dmax;
    E_Float DX, DY, DZ;
    E_Int NBIN;

    std::map<E_Int, std::set<E_Int>> bin_faces;

    Mesh();

    Mesh(const char *fname);

    Mesh(K_FLD::FldArrayI &cn, E_Float *X, E_Float *Y, E_Float *Z, E_Int npts);

    void make_skin();

    void make_bbox();

    void hash_skin();

    void make_point_faces();

    void make_edges();

    inline bool face_is_quad(E_Int face) const { return F[face].size() == 4; }
    
    inline bool face_is_tri(E_Int face) const { return F[face].size() == 3; }

    void write_ngon(const char *fname);

    void write_faces(const char *fname, const std::vector<E_Int> &faces);

    bool is_point_inside(E_Float px, E_Float py, E_Float pz);

    Mesh reconstruct_after_smesh_adaptation(const Smesh &Mf, E_Int patchc);

    // Adaptation
    void init_adaptation_data();

    std::vector<E_Int> smooth_ref_data(
        const std::map<E_Int, std::vector<E_Int>> &sensor);

    std::vector<E_Int> prepare_for_refinement(
        const std::vector<E_Int> &ref_data);

    std::set<E_Int> factive;
    std::map<E_Int, std::vector<E_Int>> fchildren;
    std::vector<E_Int> flevel;

    std::map<UEdge, E_Int> ecenter;

    size_t refine(std::set<E_Int> &mpatch, std::set<E_Int> &spatch, Mesh &S);

    std::vector<pointFace> locate(E_Float x, E_Float y,
        const std::set<E_Int> &patch) const;
    
    inline bool face_is_active(E_Int face) const
    { return factive.find(face) != factive.end(); }

    bool faces_are_dups(E_Int mface, E_Int sface, const Mesh &S);

    void refine_faces(const std::vector<E_Int> &ref_faces);

    bool face_contains_sface(E_Int face, E_Int sface, const Mesh &S) const;

    void resize_point_data(size_t nref_faces);

    void resize_face_data(size_t nref_faces);

    void refine_quad(E_Int quad);

    void refine_tri(E_Int tri);

    void refine_edge(const UEdge &edge);

    E_Int face_contains_point(E_Int face, E_Float x, E_Float y) const;

    Mesh extract_conformized();

    void get_fleaves(E_Int face, std::vector<E_Int> &fleaves);

    PyObject *export_karray();

    /* TOPO */

    E_Int orient_skin(E_Int normal_direction);

    void flag_and_get_external_faces(std::vector<E_Int> &fflags,
        std::vector<E_Int> &efaces);
    
    void extract_nface_of_kept_pgs(const std::vector<bool> &kept_pgs,
        std::vector<E_Int> &NFACE, std::vector<E_Int> &cxadj,
        std::vector<E_Int> &cells);
    
    void flag_marked_external_cells(const std::vector<E_Int> &cells,
        const std::vector<E_Int> &fflags, std::vector<E_Int> &cflags);
    
    void flag_all_external_cells(const std::vector<E_Int> &fflags,
        std::vector<E_Int> &cflags);
    
    E_Int orient_boundary(E_Int ncells, E_Int *efadj, E_Int *exadj, E_Int nefaces,
        E_Int *fneis, E_Int *efaces, std::vector<E_Int> &forient,
        const std::vector<E_Int> &cflags, const std::vector<E_Int> &fflags,
        E_Int *cells, E_Int normal_direction);
    
    void compute_cell_volume(E_Int cell, E_Float &vol, E_Int refIdx);
};

int meshes_mutual_refinement(Mesh &M, Mesh &S);