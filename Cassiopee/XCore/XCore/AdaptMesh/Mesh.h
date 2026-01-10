/*    
    Copyright 2013-2026 ONERA.

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

#include "Mpi.h"
#include <map>
#include <array>

#include "xcore.h"
#include "common/common.h"
#include "Quad.h"
#include "Hexa.h"

#define HEXA 0
#define TETRA 1
#define PENTA 2
#define PYRA 3

#define QUAD 0
#define TRI 1

#define FACE_UNTOUCHED 0
#define FACE_REFINED 1
#define FACE_NEW 2

#define ISO 0
#define DIR 1

#define DIR_ISO 0
#define DIR_X 1
#define DIR_Y 2

struct Karray;
struct SkinGraph;
struct Point;
struct Vec3f;
struct ArrayI;

struct BPatch {
    E_Int gid;
    E_Int nf;
    E_Int *pf;
    char *name;
    char *type;
};

struct PPatch {
    E_Int nei;
    E_Int nf;
    E_Int *pf;
    E_Int *pn;

    E_Int *sbuf_i;
    E_Int *rbuf_i;
};

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

    bool operator==(const UEdge &E) const
    {
        return p == E.p && q == E.q;
    }

    bool operator!=(const UEdge &E) const
    {
        return !(*this == E);
    }
};

struct Mesh {
    /* Base data */

    E_Int np;
    E_Float *X, *Y, *Z;

    E_Int nf;
    E_Int *faces;
    E_Int *frange;
    E_Int *fstride;

    E_Int ne;
    UEdge *edges; // edge points
    E_Int *fedg; // face edges

    E_Int *xedf; // edge faces offset
    E_Int *E2F; // edge faces

    E_Int nc;
    E_Int *cells; // nface equivalent
    E_Int *crange; // 1, 2 or 4 by cell face
    E_Int *cstride; // HEXA:6, TETRA:4, PENTA:5, PYRA:5

    E_Int *owner;
    E_Int *neigh;

    E_Int nconnex;

    /* Boundary */

    E_Int nbp;
    BPatch *bps;
    
    std::map<E_Int, E_Int> face_to_bpatch;

    /* Adaptation */

    E_Float *mode_2D;

    std::map<UEdge, E_Int> ecenter;

    E_Int *cref;
    E_Int *fref;
    E_Int *fpattern;

    E_Int *clevel;
    E_Int *flevel;
    E_Int *elevel;

    E_Int *ctype;
    E_Int *ftype;

    std::map<E_Int, std::vector<E_Int>> fchildren;
    std::map<E_Int, std::vector<E_Int>> cchildren;

    E_Int *fparent;

    E_Int nc_old;
    E_Int nf_old;

    E_Int gnf_old;

    /* Comm */

    E_Int npp;
    PPatch *pps;

    std::map<E_Int, E_Int> face_to_ppatch;

    /* Parallel */

    int pid;
    int npc;
    int nrq;
    MPI_Request *reqs;

    E_Int *l2gc;
    E_Int *l2gf;

    std::map<E_Int, E_Int> g2lc;
    std::map<E_Int, E_Int> g2lf;

    E_Int *xneis;
    E_Int *cneis;
    
    E_Int *ctag;
    E_Int *ftag;
    E_Int *ptag;

    Mesh();
};

/* Topo */

E_Int Mesh_set_cells_for_2D(Mesh *M);

E_Int Mesh_set_face_types(Mesh *M);

E_Int Mesh_set_cell_types(Mesh *M);

E_Int Mesh_set_orientation(Mesh *M);

E_Int Mesh_is_face_aligned_with_vec3(Mesh *M, E_Int fid, E_Float *v3);

/* Connectivity */

void Mesh_make_edge_connectivity(Mesh *M);

void Mesh_update_global_cell_ids(Mesh *M);

void Mesh_update_bpatches(Mesh *M);

void Mesh_update_ppatches(Mesh *M);

void Mesh_update_global_face_ids(Mesh *M);

E_Int Mesh_get_global_face_count(Mesh *M);

void Mesh_get_cneis(Mesh *M, E_Int cid, E_Int &nn, E_Int neis[24]);

inline
E_Int Mesh_get_cnei(Mesh *M, E_Int cid, E_Int fid)
{
    assert(cid == M->owner[fid] || cid == M->neigh[fid]);
    return (M->owner[fid] == cid) ? M->neigh[fid] : M->owner[fid];
}

inline
E_Int *Mesh_get_cell(Mesh *M, E_Int cid)
{
    return &M->cells[24*cid];
}

inline
E_Int *Mesh_get_crange(Mesh *M, E_Int cid)
{
    return &M->crange[6*cid];
}

inline
E_Int Mesh_get_sizeNFace(Mesh *M)
{
    E_Int sizeNFace = 0;

    for (E_Int i = 0; i < M->nc; i++) {
        E_Int *crange = Mesh_get_crange(M, i);
        for (E_Int j = 0; j < M->cstride[i]; j++) {
            sizeNFace += crange[j];
        }
    }
    return sizeNFace;
}

inline
E_Int *Mesh_get_face(const Mesh *M, E_Int fid)
{
    return &M->faces[8*fid];
}

inline
E_Int *Mesh_get_frange(const Mesh *M, E_Int fid)
{
    return &M->frange[4*fid];
}

inline
E_Int *Mesh_get_fedges(const Mesh *M, E_Int fid)
{
    return &M->fedg[8*fid];
}

inline
E_Int Mesh_get_sizeNGon(Mesh *M)
{
    E_Int sizeNGon = 0;

    for (E_Int i = 0; i < M->nf; i++) {
        E_Int *frange = Mesh_get_frange(M, i);
        for (E_Int j = 0; j < M->fstride[i]; j++) {
            sizeNGon += frange[j];
        }
    }
    return sizeNGon;
}


inline
void Mesh_get_fpoints(const Mesh *M, E_Int fid, E_Int &np, E_Int pts[8])
{
    np = 0;
    E_Int *face = Mesh_get_face(M, fid);
    E_Int *frange = Mesh_get_frange(M, fid);
    
    for (E_Int i = 0; i < M->fstride[fid]; i++) {
        E_Int *pn = face + 2*i;

        for (E_Int j = 0; j < frange[i]; j++) {
            pts[np++] = pn[j];
        }
    }
}

inline
void Mesh_reverse_face_points(Mesh *M, E_Int fid)
{
    E_Int *face = Mesh_get_face(M, fid);
    E_Int *frange = Mesh_get_frange(M, fid);
    E_Int fstride = M->fstride[fid];
    assert(fstride == 4);
    std::reverse(face+1, face+2*fstride);
    std::reverse(frange, frange+fstride);
}

inline E_Int Mesh_face_is_iface(Mesh *M, E_Int lfid)
{
    return M->owner[lfid] != -1 && M->neigh[lfid] != -1;
}

inline E_Int Mesh_face_is_pface(Mesh *M, E_Int lfid)
{
    return M->face_to_ppatch.find(lfid) != M->face_to_ppatch.end();
}

inline E_Int Mesh_face_is_bface(Mesh *M, E_Int lfid)
{
    return M->face_to_bpatch.find(lfid) != M->face_to_bpatch.end();
}

E_Int Mesh_build_pe(Mesh *M);

void Mesh_make_cell_cells(Mesh *M);

inline E_Int Mesh_get_reorient(Mesh *M, E_Int face, E_Int cell, E_Int normalIn)
{
    assert(M->owner[face] == cell || M->neigh[face] == cell);
    E_Int ret = (M->neigh[face] == cell && normalIn == 0) ||
              (M->owner[face] == cell && normalIn == 1);
    return ret;
}


/* Clean-up */

void Mesh_reset_base_data(Mesh *M);

void Mesh_reset_boundary_data(Mesh *M);

void Mesh_reset_comm_data(Mesh *M);

void Mesh_reset_adaptation_data(Mesh *M);

void Mesh_reset_parallel_data(Mesh *M);

void Mesh_reset_tags(Mesh *M);

void Mesh_free(Mesh *M);


/* IO */

Mesh *Mesh_from_Karray(Karray *karray);

PyObject *Mesh_export_karray(Mesh *M, E_Int conformize);


/* Parallel */

inline
void Mesh_comm_waitall(Mesh *M)
{
    MPI_Waitall(M->nrq, M->reqs, MPI_STATUSES_IGNORE);
    M->nrq = 0;
}

E_Int Mesh_load_balance(Mesh *M);

/* Adaptation */

void Mesh_smooth_skin_ref_data(Mesh *M, const SkinGraph *skin_graph,
    E_Int *fdat);

void Mesh_smooth_cell_refinement_data(Mesh *M);

E_Int Mesh_smooth_cref(Mesh *M);

void Mesh_get_ref_entities(Mesh *M, std::vector<E_Int> &ref_cells,
    std::vector<E_Int> &ref_faces, std::set<UEdge> &ref_edges);

void Mesh_resize_for_refinement(Mesh *M, const std::vector<E_Int> &ref_cells,
    const std::vector<E_Int> &ref_faces, const std::set<UEdge> &ref_edges);

void Mesh_resize_face_data(Mesh *M, E_Int new_nf);

void Mesh_resize_cell_data(Mesh *M, E_Int new_nc);

void Mesh_sort_ref_entities_by_level(Mesh *M,
    std::vector<E_Int> &ref_cells, std::vector<E_Int> &ref_faces,
    std::set<UEdge> &ref_edges);

void Mesh_refine(Mesh *M, std::vector<E_Int> &ref_cells,
    std::vector<E_Int> &ref_faces, std::set<UEdge> &ref_edges);

E_Int Mesh_conformize_cell_face(Mesh *M, E_Int cid, E_Int fid, E_Int fpos, E_Int nf);

void Mesh_conformize_face_edge(Mesh *M);

void Mesh_refine_iso(Mesh *M, std::vector<E_Int> &ref_cells,
    std::vector<E_Int> &ref_faces, std::set<UEdge> &ref_edges);

void Mesh_refine_dir(Mesh *M, std::vector<E_Int> &ref_cells,
    std::vector<E_Int> &ref_faces, std::set<UEdge> &ref_edges);

inline
void update_shell_pe(E_Int parent, Mesh *M)
{
    const auto &children = M->cchildren.at(parent);

    for (E_Int cid : children) {
        E_Int *child = Mesh_get_cell(M, cid);

        for (E_Int j = 0; j < M->cstride[cid]; j++) {
            E_Int face = child[4*j];
            
            if      (M->owner[face] == parent) M->owner[face] = cid;
            else if (M->neigh[face] == parent) M->neigh[face] = cid;
        }
    }
}

inline
void Mesh_update_face_range_and_stride(Mesh *M, E_Int quad, E_Int fpos, E_Int nchild)
{
    E_Int *frange = Mesh_get_frange(M, quad);
    for (E_Int i = 0; i < M->fstride[quad]; i++) {
        frange[i] = 1;
    }

    for (E_Int i = 0; i < nchild; i++) {
        E_Int child = fpos + i;

        M->fstride[child] = M->fstride[quad];

        frange = Mesh_get_frange(M, child);
        for (E_Int j = 0; j < M->fstride[child]; j++) {
            frange[j] = 1;
        }
    }
}

struct point {
    E_Float x, y, z;
};

void write_points(const char *fname, const std::vector<point> &points);

void write_point(const char *fname, E_Float x, E_Float y, E_Float z);

inline
void Mesh_refine_or_get_edge_center(Mesh *M, E_Int p, E_Int q, E_Int &node)
{
    UEdge E(p, q);
    auto it = M->ecenter.find(E);
    if (it == M->ecenter.end()) {
        M->X[M->np] = 0.5 * (M->X[p] + M->X[q]);
        M->Y[M->np] = 0.5 * (M->Y[p] + M->Y[q]);
        M->Z[M->np] = 0.5 * (M->Z[p] + M->Z[q]);
        /*if (M->np == 29) {
            write_point("p", M->X[p], M->Y[p], M->Z[p]);
            write_point("q", M->X[q], M->Y[q], M->Z[q]);
            write_point("c", M->X[29], M->Y[29], M->Z[29]);
            write_point("o", M->X[9], M->Y[9], M->Z[9]);
        }*/
        node = M->np;
        M->ecenter[E] = M->np;
        M->np++;
    } else {
        node = it->second;
    }
}

inline
void refine_cell_dir(E_Int cell, Mesh *M)
{
    switch (M->ctype[cell]) {
        case HEXA:
            H18_refine(cell, M);
            break;
        default:
            assert(0);
            break;
    }
}

inline
void refine_face_dir(E_Int face, E_Int pattern, Mesh *M)
{
    switch (M->ftype[face]) {
        case QUAD: {
            if (pattern == ISO) Q9_refine(face, M);
            else Q6_refine(face, M);
            break;
        }
        default:
            assert(0);
            break;
    }
}

inline
E_Int get_face_pattern(E_Int fid, Mesh *M)
{
    E_Int own = M->owner[fid];
    E_Int fpos = -1;
    E_Int *cell = Mesh_get_cell(M, own);
    E_Int *crange = Mesh_get_crange(M, own);

    for (E_Int j = 0; j < M->cstride[own] && fpos == -1; j++) {
        E_Int *pf = cell + 4*j;
        for (E_Int k = 0; k < crange[j]; k++) {
            E_Int face = pf[k];
            if (face == fid) {
                fpos = j;
                assert(k == 0);
                break;
            }
        }
    }

    assert(fpos != -1);

    if (fpos == 0 || fpos == 1) return ISO;
    return DIR;
}

void Mesh_triangulate_face(Mesh *M, E_Int fid);
void Mesh_triangulate_faces(Mesh *M, E_Int *faces, E_Int nf);

void Mesh_face_to_prism(Mesh *M, E_Int fid);
void Mesh_generate_prisms(Mesh *M, E_Int *faces, E_Int nf);

/* Extract */

void Mesh_extract_skin(const Mesh *M, E_Int *count, E_Int **skin);

void Mesh_make_skin_connectivity(const Mesh *M, SkinGraph *skin_graph);

void Mesh_make_skin_graph(const Mesh *M, SkinGraph *skin_graph);

void Mesh_make_face_centers(const Mesh *M, const E_Int nf, const E_Int *skin,
    Vec3f *fc);

void Mesh_extract_points_from_ftag(const Mesh *M, ArrayI *pids);

void Mesh_SkinGraph_compute_normals(const Mesh *M, SkinGraph *skin_graph);

/* Locate */

bool Mesh_point_in_tri(const Mesh *M, const Point *p, E_Int tid);

bool Mesh_point_in_quad(const Mesh *M, const Point *p, E_Int qid);

bool Mesh_point_in_face(const Mesh *M, const Point *p, E_Int fid);