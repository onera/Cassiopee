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

#include <mpi.h>
#include <map>
#include <array>

#include "xcore.h"
#include "../common/common.h"

#define HEXA 0
#define TETRA 1
#define PENTA 2
#define PYRA 3

#define QUAD 0
#define TRI 1

#define FACE_UNTOUCHED 0
#define FACE_REFINED 1
#define FACE_NEW 2

struct Karray;

struct BPatch {
    Int gid;
    Int nf;
    Int *pf;
    char *name;
    char *type;
};

struct PPatch {
    Int nei;
    Int nf;
    Int *pf;
    Int *pn;

    Int *sbuf_i;
    Int *rbuf_i;
};

struct UEdge {
    Int p, q;

    UEdge(Int P, Int Q)
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

    Int np;
    Float *X, *Y, *Z;

    Int nf;
    Int *faces;
    Int *frange;
    Int *fstride;

    Int ne;
    UEdge *edges; // edge points
    Int *fedg; // face edges

    Int *xedf; // edge faces offset
    Int *E2F; // edge faces

    Int nc;
    Int *cells; // nface equivalent
    Int *crange; // 1, 2 or 4 by cell face
    Int *cstride; // HEXA:6, TETRA:4, PENTA:5, PYRA:5

    Int *owner;
    Int *neigh;

    Int nconnex;

    /* Boundary */

    Int nbp;
    BPatch *bps;
    
    std::map<Int, Int> face_to_bpatch;

    /* Adaptation */

    Float *mode_2D;

    std::map<UEdge, Int> ecenter;

    Int *cref;
    Int *fref;

    Int *clevel;
    Int *flevel;
    Int *elevel;

    Int *ctype;
    Int *ftype;

    std::map<Int, std::vector<Int>> fchildren;
    std::map<Int, std::vector<Int>> cchildren;

    Int *fparent;

    Int nc_old;
    Int nf_old;

    Int gnf_old;

    /* Comm */

    Int npp;
    PPatch *pps;

    std::map<Int, Int> face_to_ppatch;

    /* Parallel */

    Int pid;
    Int npc;
    Int nrq;
    MPI_Request *reqs;

    Int *l2gc;
    Int *l2gf;

    std::map<Int, Int> g2lc;
    std::map<Int, Int> g2lf;

    Int *xneis;
    Int *cneis;

    Mesh();
};

/* Topo */

Int Mesh_set_cells_for_2D(Mesh *M);

void Mesh_get_cneis(Mesh *M, Int cid, Int &nn, Int neis[24]);

Int Mesh_set_face_types(Mesh *M);

Int Mesh_set_cell_types(Mesh *M);

Int Mesh_set_orientation(Mesh *M);

Int Mesh_is_face_aligned_with_vec3(Mesh *M, Int fid, Float *v3);

/* Connectivity */

void Mesh_make_edge_connectivity(Mesh *M);

void Mesh_update_global_cell_ids(Mesh *M);

void Mesh_update_bpatches(Mesh *M);

void Mesh_update_ppatches(Mesh *M);

void Mesh_update_global_face_ids(Mesh *M);

inline
Int *Mesh_get_cell(Mesh *M, int cid)
{
    return &M->cells[24*cid];
}

inline
Int *Mesh_get_crange(Mesh *M, int cid)
{
    return &M->crange[6*cid];
}

inline
Int Mesh_get_sizeNFace(Mesh *M)
{
    Int sizeNFace = 0;

    for (Int i = 0; i < M->nc; i++) {
        Int *crange = Mesh_get_crange(M, i);
        for (Int j = 0; j < M->cstride[i]; j++) {
            sizeNFace += crange[j];
        }
    }
    return sizeNFace;
}

inline
Int *Mesh_get_face(Mesh *M, int fid)
{
    return &M->faces[8*fid];
}

inline
Int *Mesh_get_frange(Mesh *M, int fid)
{
    return &M->frange[4*fid];
}

inline
Int *Mesh_get_fedges(Mesh *M, Int fid)
{
    return &M->fedg[8*fid];
}

inline
Int Mesh_get_sizeNGon(Mesh *M)
{
    Int sizeNGon = 0;

    for (Int i = 0; i < M->nf; i++) {
        Int *frange = Mesh_get_frange(M, i);
        for (Int j = 0; j < M->fstride[i]; j++) {
            sizeNGon += frange[j];
        }
    }
    return sizeNGon;
}


inline
void Mesh_get_fpoints(Mesh *M, Int fid, Int &np, Int pts[8])
{
    np = 0;
    Int *face = Mesh_get_face(M, fid);
    Int *frange = Mesh_get_frange(M, fid);
    
    for (Int i = 0; i < M->fstride[fid]; i++) {
        Int *pn = face + 2*i;

        for (Int j = 0; j < frange[i]; j++) {
            pts[np++] = pn[j];
        }
    }
}

inline
void Mesh_reverse_face_points(Mesh *M, Int fid)
{
    Int *face = Mesh_get_face(M, fid);
    Int *frange = Mesh_get_frange(M, fid);
    Int fstride = M->fstride[fid];
    assert(fstride == 4);
    std::reverse(face+1, face+2*fstride);
    std::reverse(frange, frange+fstride);
}

inline Int Mesh_face_is_iface(Mesh *M, Int lfid)
{
    return M->owner[lfid] != -1 && M->neigh[lfid] != -1;
}

inline Int Mesh_face_is_pface(Mesh *M, Int lfid)
{
    return M->face_to_ppatch.find(lfid) != M->face_to_ppatch.end();
}

inline Int Mesh_face_is_bface(Mesh *M, Int lfid)
{
    return M->face_to_bpatch.find(lfid) != M->face_to_bpatch.end();
}

void Mesh_make_cell_cells(Mesh *M);

inline Int Mesh_get_reorient(Mesh *M, Int face, Int cell, Int normalIn)
{
    assert(M->owner[face] == cell || M->neigh[face] == cell);
    Int ret = (M->neigh[face] == cell && normalIn == 0) ||
              (M->owner[face] == cell && normalIn == 1);
    return ret;
}


/* Clean-up */

void Mesh_reset_base_data(Mesh *M);

void Mesh_reset_boundary_data(Mesh *M);

void Mesh_reset_comm_data(Mesh *M);

void Mesh_reset_adaptation_data(Mesh *M);

void Mesh_reset_parallel_data(Mesh *M);

void Mesh_free(Mesh *M);


/* IO */

Mesh *Mesh_from_Karray(Karray *karray);

PyObject *Mesh_export_karray(Mesh *M, int conformize);


/* Parallel */

void Mesh_comm_waitall(Mesh *M);

Int Mesh_load_balance(Mesh *M);


/* Adaptation */

Int Mesh_smooth_cref(Mesh *M);

void Mesh_get_ref_entities(Mesh *M, std::vector<Int> &ref_cells,
    std::vector<Int> &ref_faces, std::set<UEdge> &ref_edges);

void Mesh_resize_for_refinement(Mesh *M, const std::vector<Int> &ref_cells,
    const std::vector<Int> &ref_faces, const std::set<UEdge> &ref_edges);

void Mesh_sort_ref_entities_by_level(Mesh *M,
    std::vector<Int> &ref_cells, std::vector<Int> &ref_faces,
    std::set<UEdge> &ref_edges);

void Mesh_refine(Mesh *M, std::vector<Int> &ref_cells,
    std::vector<Int> &ref_faces, std::set<UEdge> &ref_edges);

int Mesh_conformize_cell_face(Mesh *M, Int cid, Int fid, Int fpos, Int nf);

void Mesh_conformize_face_edge(Mesh *M);

void Mesh_refine_iso(Mesh *M, std::vector<Int> &ref_cells,
    std::vector<Int> &ref_faces, std::set<UEdge> &ref_edges);

void Mesh_refine_dir(Mesh *M, std::vector<Int> &ref_cells,
    std::vector<Int> &ref_faces, std::set<UEdge> &ref_edges);

inline
void Mesh_update_face_range_and_stride(Mesh *M, Int quad, Int fpos, Int nchild)
{
    Int *frange = Mesh_get_frange(M, quad);
    for (Int i = 0; i < M->fstride[quad]; i++) {
        frange[i] = 1;
    }

    for (Int i = 0; i < nchild; i++) {
        Int child = fpos + i;

        M->fstride[child] = M->fstride[quad];

        frange = Mesh_get_frange(M, child);
        for (Int j = 0; j < M->fstride[child]; j++) {
            frange[j] = 1;
        }
    }
}

inline
void Mesh_refine_or_get_edge_center(Mesh *M, Int p, Int q, Int &node)
{
    UEdge E(p, q);
    auto it = M->ecenter.find(E);
    if (it == M->ecenter.end()) {
        M->X[M->np] = 0.5 * (M->X[p] + M->X[q]);
        M->Y[M->np] = 0.5 * (M->Y[p] + M->Y[q]);
        M->Z[M->np] = 0.5 * (M->Z[p] + M->Z[q]);
        node = M->np;
        M->ecenter[E] = M->np;
        M->np++;
    } else {
        node = it->second;
    }
}