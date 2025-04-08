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
#include <list>
#include <unordered_map>

#include "point.h"
#include "xcore.h"
#include "common/common.h"
#include "triangle.h"
#include "AABB.h"
#include "smesh.h"

#define OUT 0
#define IN 1

struct Ray;
struct Karray;

struct UEdge {
    E_Int p, q;

    UEdge(E_Int P, E_Int Q)
    {
        p = std::min(P, Q);
        q = std::max(P, Q);
    }

    bool operator<(const UEdge &E) const
    {
        return (p < E.p) || (p == E.p && q < E.q);
    }
};

struct Sgraph {
    std::vector<E_Int> xadj;
    std::vector<E_Int> fpts;
    std::vector<E_Int> fadj;
};

struct Py_BC {
    E_Int size;
    E_Int *ptr;
};

struct IMesh {
    E_Int np, ne, nf, nc;

    std::vector<E_Float> X, Y, Z;
    std::vector<std::vector<E_Int>> F;

    std::vector<std::vector<E_Int>> P2F;
    std::vector<std::vector<E_Int>> P2E;

    std::vector<o_edge> E;
    std::vector<std::vector<E_Int>> E2F;
    std::vector<std::vector<E_Int>> F2E;

    std::vector<std::vector<E_Int>> C;

    std::vector<E_Int> skin;
    std::vector<E_Int> owner;
    std::vector<E_Int> neigh;

    E_Float NEAR_VERTEX_TOL = 1e-3;
    E_Float NEAR_EDGE_TOL = 1e-3;

    E_Float xmin, ymin, zmin;
    E_Float xmax, ymax, zmax;
    E_Int NX, NY, NZ;
    E_Int NXY, NXYZ;
    E_Float HX, HY, HZ;

    std::vector<std::vector<E_Int>> bin_faces;

    std::set<E_Int> patch;
    std::vector<E_Int> ftag;
    std::vector<E_Float> ptag;
    std::vector<E_Float> ctag;

    Smesh Mf;

    E_Float get_min_edge_length() const;
    
    void set_tolerances(E_Float near_vertex_tol, E_Float near_edge_tol)
    {
        NEAR_VERTEX_TOL = near_vertex_tol;
        NEAR_EDGE_TOL = near_edge_tol;
    }

    inline E_Int get_voxel(E_Int I, E_Int J, E_Int K) const
    {
        return I + J * NX + (NX * NY) * K;
    }

    AABB AABB_face(const std::vector<E_Int> &pn) const;
    
    void triangulate_face_set(bool propagate = true);

    void triangulate(const std::vector<E_Int> &fids);

    void triangulate_skin(std::unordered_map<E_Int, E_Int> &fid_to_bc);

    void triangulate_skin();

    size_t refine_slave(const IMesh &master);
    
    void hash_patch();

    IMesh();

    IMesh(const Karray &karray);

    //IMesh(K_FLD::FldArrayI &cn, E_Float *X, E_Float *Y, E_Float *Z, E_Int npts);

    void make_patch(E_Int *faces, E_Int nfaces);

    void make_skin();

    void make_bbox();

    void hash_skin();

    void make_point_faces();

    void make_edges();

    void extract_edge_points(E_Int a, E_Int b, std::list<E_Int> &pts);

    inline bool face_is_quad(E_Int face) const { return F[face].size() == 4; }
    
    inline bool face_is_tri(E_Int face) const { return F[face].size() == 3; }

    void write_ngon(const char *fname);

    void write_ngon(const char *fname, const std::vector<E_Int> &faces) const;

    void write_ngon(const char *fname, const std::set<E_Int> &fset) const;

    void write_face(const char *fname, E_Int fid) const;

    bool is_point_inside(E_Float px, E_Float py, E_Float pz) const;

    IMesh reconstruct_after_smesh_adaptation(const Smesh &Mf, E_Int patchc);

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

    std::set<E_Int> faces_to_tri;

    size_t refine(const IMesh &S);

    inline bool face_is_active(E_Int face) const
    { return factive.find(face) != factive.end(); }

    bool faces_are_dups(E_Int mface, E_Int sface, const IMesh &S);

    void refine_faces(const std::vector<E_Int> &ref_faces);

    bool face_contains_sface(E_Int face, E_Int sface, const IMesh &S) const;

    void resize_point_data(size_t nref_faces);

    void resize_face_data(size_t nref_faces);

    void refine_quad(E_Int quad);

    void refine_tri(E_Int tri);

    void refine_edge(const UEdge &edge);

    E_Int face_contains_point(E_Int face, E_Float x, E_Float y, E_Float z) const;

    IMesh extract_conformized();

    void get_fleaves(E_Int face, std::vector<E_Int> &fleaves);

    PyObject *export_karray(E_Int remove_periodic = 0) const;

    PyObject *export_karray_orig() const;

    PyObject *export_karray_periodic() const;

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

E_Int meshes_mutual_refinement(IMesh &M, IMesh &S);
