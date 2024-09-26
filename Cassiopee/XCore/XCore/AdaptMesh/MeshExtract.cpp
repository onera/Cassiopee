#include "Mesh.h"
#include "Skin.h"
#include "common/mem.h"
#include "Array.h"
#include "FaceSort.h"
#include "Vec.h"

void Mesh_extract_points_from_ftag(const Mesh *M, ArrayI *pids)
{
    E_Int *ptag = (E_Int *)XMALLOC(M->np * sizeof(E_Int));
    memset(ptag, 0, M->np * sizeof(E_Int));
    pids->count = 0;

    for (E_Int fid = 0; fid < M->nf; fid++) {
        if (M->ftag[fid] != 1) continue;
        E_Int *face = Mesh_get_face(M, fid);
        E_Int *frange = Mesh_get_frange(M, fid);
        for (E_Int i = 0; i < M->fstride[fid]; i++) {
            E_Int *pn = face + 2*i;
            for (E_Int j = 0; j < frange[i]; j++) {
                E_Int pid = pn[j];
                pids->count += (ptag[pid] == 0);
                ptag[pid] = 1;
            }
        }
    }

    pids->ptr = (E_Int *)XMALLOC(pids->count * sizeof(E_Int));
    E_Int *ptr = pids->ptr;

    for (E_Int pid = 0; pid < M->np; pid++) {
        if (ptag[pid] == 1)
            *ptr++ = pid;
    }

    XFREE(ptag);
}

void Mesh_extract_faces_from_ftag(const Mesh *M, FaceSort **mfaces,
    E_Int *mcount)
{
    E_Int count = 0;

    for (E_Int i = 0; i < M->nf; i++) {
        count += (M->ftag[i] == 1);
    }

    *mfaces = (FaceSort *)XMALLOC(count * sizeof(FaceSort));
    
    *mcount = count;

    count = 0;

    for (E_Int i = 0; i < M->nf; i++) {
        if (M->ftag[i] == 1) {
            assert(M->ftype[i] == TRI);
            (*mfaces)[count++].fid = i;
        }
    }

    assert(*mcount == count);
}

void Mesh_extract_skin(const Mesh *M, E_Int *count, E_Int **skin)
{
    *count = 0;
    
    for (E_Int fid = 0; fid < M->nf; fid++) {
        *count += (M->neigh[fid] == -1);
    }

    *skin = (E_Int *)XMALLOC(*count * sizeof(E_Int));
    E_Int *ptr = *skin;

    for (E_Int fid = 0; fid < M->nf; fid++) {
        if (M->neigh[fid] == -1)
            *ptr++ = fid;
    }
}

void Mesh_make_skin_connectivity(const Mesh *M, SkinGraph *skin_graph)
{
    // Count
    skin_graph->xadj = (E_Int *)XMALLOC((skin_graph->nf+1) * sizeof(E_Int));
    E_Int *xadj = skin_graph->xadj;
    xadj[0] = 0;

    for (E_Int i = 0; i < skin_graph->nf; i++) {
        E_Int fid = skin_graph->skin[i];
        const E_Int *frange = Mesh_get_frange(M, fid);
        xadj[i+1] = 0;
        for (E_Int j = 0; j < M->fstride[fid]; j++)
            xadj[i+1] += frange[j];
        xadj[i+1] += xadj[i];
    }

    skin_graph->fpts = (E_Int *)XMALLOC(xadj[skin_graph->nf] * sizeof(E_Int)); 

    // Populate
    E_Int *ptr = skin_graph->fpts;

    for (E_Int i = 0; i < skin_graph->nf; i++) {
        E_Int fid = skin_graph->skin[i];
        const E_Int *face = Mesh_get_face(M, fid);
        const E_Int *frange = Mesh_get_frange(M, fid);
        for (E_Int j = 0; j < M->fstride[fid]; j++) {
            const E_Int *pn = face + 2*j;
            for (E_Int k = 0; k < frange[j]; k++)
                *ptr++ = pn[k];
        }
    }
}

void Mesh_make_skin_graph(const Mesh *M, SkinGraph *skin_graph)
{
    Mesh_extract_skin(M, &skin_graph->nf, &skin_graph->skin);
    Mesh_make_skin_connectivity(M, skin_graph);
    SkinGraph_make_skin_neighbours(skin_graph);
}

void Mesh_make_face_centers(const Mesh *M, const E_Int nf, const E_Int *skin,
    Vec3f *fc)
{
    for (E_Int i = 0; i < nf; i++) {
        E_Int fid = skin[i];
        const E_Int *face = Mesh_get_face(M, fid);
        const E_Int *frange = Mesh_get_frange(M, fid);
        fc[i].x = fc[i].y = fc[i].z = 0.0;
        E_Int np = 0;
        for (E_Int j = 0; j < M->fstride[fid]; j++) {
            const E_Int *pn = face + 2*j;
            for (E_Int k = 0; k < frange[j]; k++) {
                fc[i].x += M->X[pn[k]];
                fc[i].y += M->Y[pn[k]];
                fc[i].z += M->Z[pn[k]];
                np++;
            }
        }
        fc[i].x /= np; fc[i].y /= np; fc[i].z /= np;
    }
}

void Mesh_SkinGraph_compute_normals(const Mesh *M, SkinGraph *skin_graph)
{
    assert(!skin_graph->fnml);
    skin_graph->fnml = (Vec3f *)XMALLOC(skin_graph->nf * sizeof(Vec3f));

    for (E_Int i = 0; i < skin_graph->nf; i++) {
        E_Int fid = skin_graph->skin[i];
        E_Int *face = Mesh_get_face(M, fid);
        E_Int A = face[0], B = face[2], C = face[4];
        E_Float U[3] = {M->X[B]-M->X[A], M->Y[B]-M->Y[A], M->Z[B]-M->Z[A]};
        E_Float V[3] = {M->X[C]-M->X[A], M->Y[C]-M->Y[A], M->Z[C]-M->Z[A]};
        Vec3f *N = &skin_graph->fnml[i];
        K_MATH::cross(U, V, (E_Float *)N);
        E_Float NORM = K_MATH::norm((E_Float *)N, 3);
        N->x /= NORM;
        N->y /= NORM;
        N->z /= NORM;
    }
}