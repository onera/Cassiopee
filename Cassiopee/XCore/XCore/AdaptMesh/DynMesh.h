#pragma once

#include <vector>
#include <map>
#include <set>

#include "xcore.h"
#include "common/common.h"

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

    std::vector<E_Int> skin;

    std::vector<E_Int> ftag;

    std::vector<E_Int> fref;

    DynMesh();

    DynMesh(Karray *karray);

    void extract_points_from_ftag(ArrayI *points);

    void make_skin_graph(SkinGraph *skin_graph);

    void make_face_centers(const E_Int NF, const E_Int *skin,
        Vec3f *fc);

    void smooth_skin_ref_data(SkinGraph *skin_graph, E_Int *fdat);

    bool point_in_tri(const Point *p, E_Int tid) const;

    bool point_in_quad(const Point *p, E_Int qid) const;

    bool point_in_face(const Point *p, E_Int fid) const;

    E_Int orient_skin(E_Int normal_direction);

    void extract_skin(E_Int *count, E_Int **skin);

    void make_skin_connectivity(SkinGraph *skin_graph);

    void make_skin_neighbours(SkinGraph *skin_graph);


    inline bool face_is_quad(E_Int face) const { return F[face].size() == 4; }
    
    inline bool face_is_tri(E_Int face) const { return F[face].size() == 3; }

    void write_ngon(const char *fname);

    void write_faces(const char *fname, const std::vector<E_Int> &faces) const;

    void write_face(const char *fname, E_Int fid) const;

    // Adaptation
    void init_adaptation_data();

    std::vector<E_Int> smooth_ref_data(
        const std::map<E_Int, std::vector<E_Int>> &sensor);

    std::vector<E_Int> prepare_for_refinement(
        const std::vector<E_Int> &ref_data);

    std::set<E_Int> factive;
    std::map<E_Int, std::vector<E_Int>> fchildren;
    std::vector<E_Int> flevel;

    std::map<uEdge, E_Int> ecenter;

    size_t refine(const DynMesh &S);
    
    inline bool face_is_active(E_Int face) const
    { return factive.find(face) != factive.end(); }

    void refine_faces(const std::vector<E_Int> &ref_faces);

    void resize_point_data(size_t nref_faces);

    void resize_face_data(size_t nref_faces);

    void refine_quad(E_Int quad);

    void refine_tri(E_Int tri);

    DynMesh extract_conformized();

    void get_fleaves(E_Int face, std::vector<E_Int> &fleaves);

    PyObject *export_karray(E_Int remove_periodic = 0);

    PyObject *export_karray_orig();

    PyObject *export_karray_periodic();

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
