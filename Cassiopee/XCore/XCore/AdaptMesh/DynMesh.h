#pragma once

#include <vector>
#include <map>
#include <set>
#include <list>

#include "xcore.h"
#include "common/common.h"
#include "TriGraph.h"

struct Karray;
struct ArrayI;
struct SkinGraph;
struct Vec3f;
struct Point;

#define OUT 0
#define IN 1

struct uEdge {
    E_Int p, q;

    uEdge(E_Int P, E_Int Q)
    {
        p = std::min(P, Q);
        q = std::max(P, Q);
    }

    bool operator<(const uEdge &E) const
    {
        return (p < E.p) || (p == E.p && q < E.q);
    }
};

struct DynMesh {
    E_Int np, ne, nf, nc;

    std::vector<E_Float> X, Y, Z;
    
    std::vector<std::vector<E_Int>> F;

    std::vector<std::vector<E_Int>> C;

    std::vector<E_Int> owner, neigh;

    std::vector<E_Int> ftag;

    std::map<E_Int, std::vector<E_Int>> fchildren;

    std::map<uEdge, E_Int> ecenter;

    TriGraph tri_graph;

    DynMesh();

    DynMesh(Karray *karray);

    void triangulate(const E_Int *faces, E_Int fcount);

    E_Int get_edge_index(E_Int nei, E_Int tri);

    void extract_points_from_ftag(ArrayI *points);

    void make_tri_graph();

    void make_face_centers(const E_Int NF, const E_Int *skin,
        Vec3f *fc);

    bool point_in_tri(const Point *p, E_Int tid) const;

    bool point_in_quad(const Point *p, E_Int qid) const;

    bool point_in_face(const Point *p, E_Int fid) const;

    E_Int orient_skin(E_Int normal_direction);

    void extract_skin(E_Int *count, E_Int **skin);

    void make_skin_connectivity(SkinGraph *skin_graph);

    void make_skin_neighbours(SkinGraph *skin_graph);

    void prepare_for_refinement(ArrayI *ref_faces);

    void refine_faces(ArrayI *ref_faces);

    inline bool face_is_tri(E_Int face) const { return F[face].size() == 3; }

    DynMesh extract_conformized();

    void write_ngon(const char *fname);

    void write_faces(const char *fname, const std::vector<E_Int> &faces) const;

    void write_face(const char *fname, E_Int fid) const;
    
    void resize_point_data(size_t nref_faces);

    void resize_face_data(size_t nref_faces);

    void refine_tri(E_Int tri);

    void get_fleaves(E_Int face, std::vector<E_Int> &fleaves);

    PyObject *export_karray();

    void extract_edge_points(E_Int a, E_Int b, std::list<E_Int> &points);

    /* TOPO */

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
