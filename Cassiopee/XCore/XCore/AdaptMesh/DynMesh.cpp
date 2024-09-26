#include "DynMesh.h"
#include "Karray.h"
#include "Array.h"
#include "Skin.h"
#include "Vec.h"
#include "Point.h"
#include "common/mem.h"

DynMesh::DynMesh(Karray *karray)
{
    np = karray->npoints();
    nf = karray->nfaces();
    nc = karray->ncells();

    X.resize(np);
    Y.resize(np);
    Z.resize(np);
    memcpy(X.data(), karray->X(), np * sizeof(E_Float));
    memcpy(Y.data(), karray->Y(), np * sizeof(E_Float));
    memcpy(Z.data(), karray->Z(), np * sizeof(E_Float));

    F.reserve(nf);

    for (E_Int i = 0; i < nf; i++) {
        E_Int np = -1;
        E_Int *pn = karray->get_face(i, np);
        std::vector<E_Int> points(np);
        for (E_Int j = 0; j < np; j++)
            points[j] = pn[j] - 1;
        F.push_back(points);
    }

    C.reserve(nc);
    for (E_Int i = 0; i < nc; i++) {
        E_Int nf = -1;
        E_Int *pf = karray->get_cell(i, nf);
        std::vector<E_Int> faces(nf);
        for (E_Int j = 0; j < nf; j++)
            faces[j] = pf[j] - 1;
        C.push_back(faces);
    }

    ftag.resize(nf, 0);

    owner.resize(nf, -1);
    neigh.resize(nf, -1);
    for (E_Int cid = 0; cid < nc; cid++) {
        const auto &pf = C[cid];
        for (E_Int fid : pf) {
            if (owner[fid] == -1) owner[fid] = cid;
            else neigh[fid] = cid;
        }
    }
}

void DynMesh::extract_points_from_ftag(ArrayI *pids)
{
    E_Int *ptag = (E_Int *)XMALLOC(np * sizeof(E_Int));
    memset(ptag, 0, np * sizeof(E_Int));
    pids->count = 0;

    for (E_Int fid = 0; fid < nf; fid++) {
        if (ftag[fid] != 1) continue;
        const auto &pn = F[fid];
        
        for (E_Int pid : pn) {
            pids->count += (ptag[pid] == 0);
            ptag[pid] = 1;
        }
    }

    pids->ptr = (E_Int *)XMALLOC(pids->count * sizeof(E_Int));
    E_Int *ptr = pids->ptr;

    for (E_Int pid = 0; pid < np; pid++) {
        if (ptag[pid] == 1)
            *ptr++ = pid;
    }

    XFREE(ptag);
}

void DynMesh::extract_skin(E_Int *count, E_Int **skin)
{
    *count = 0;
    
    for (E_Int fid = 0; fid < nf; fid++) {
        *count += (neigh[fid] == -1);
    }

    *skin = (E_Int *)XMALLOC(*count * sizeof(E_Int));
    E_Int *ptr = *skin;

    for (E_Int fid = 0; fid < nf; fid++) {
        if (neigh[fid] == -1)
            *ptr++ = fid;
    }
}

void DynMesh::make_skin_connectivity(SkinGraph *skin_graph)
{
    // Count
    skin_graph->xadj = (E_Int *)XMALLOC((skin_graph->nf+1) * sizeof(E_Int));
    E_Int *xadj = skin_graph->xadj;
    xadj[0] = 0;

    for (E_Int i = 0; i < skin_graph->nf; i++) {
        E_Int fid = skin_graph->skin[i];
        xadj[i+1] = F[fid].size();
        xadj[i+1] += xadj[i];
    }

    skin_graph->fpts = (E_Int *)XMALLOC(xadj[skin_graph->nf] * sizeof(E_Int)); 

    // Populate
    E_Int *ptr = skin_graph->fpts;

    for (E_Int i = 0; i < skin_graph->nf; i++) {
        E_Int fid = skin_graph->skin[i];
        const auto &pn = F[fid];
        for (E_Int p : pn) {
            *ptr++ = p;
        }
    }
}

void DynMesh::make_skin_graph(SkinGraph *skin_graph)
{
    extract_skin(&skin_graph->nf, &skin_graph->skin);
    make_skin_connectivity(skin_graph);
    SkinGraph_make_skin_neighbours(skin_graph);
}

void DynMesh::make_face_centers(const E_Int NF, const E_Int *skin,
    Vec3f *fc)
{
    for (E_Int i = 0; i < NF; i++) {
        fc[i].x = fc[i].y = fc[i].z = 0.0;
        E_Int fid = skin[i];
        const auto &pn = F[fid];
        for (E_Int p : pn) {
            fc[i].x += X[p];
            fc[i].y += Y[p];
            fc[i].z += Z[p];
        }
        fc[i].x /= pn.size(); fc[i].y /= pn.size(); fc[i].z /= pn.size();
    }
}

bool DynMesh::point_in_tri(const Point *p, E_Int tid) const
{
    const auto &pn = F[tid];
    E_Int A = pn[0], B = pn[1], C = pn[2];
    return Point_in_tri(p->x, p->y, p->z,
                        X[A], Y[A], Z[A],
                        X[B], Y[B], Z[B],
                        X[C], Y[C], Z[C]);
}

bool DynMesh::point_in_quad(const Point *p, E_Int qid) const
{
    // TODO(Imad): maybe compute face centers once in pre-pass
    // Star the quad into 4 triangles
    E_Float O[3] = {0.0, 0.0, 0.0};
    const auto &pn = F[qid];
    E_Int A = pn[0], B = pn[1], C = pn[2], D = pn[3];
    O[0] = (X[A] + X[B] + X[C] + X[D]) * 0.25;
    O[1] = (Y[A] + Y[B] + Y[C] + Y[D]) * 0.25;
    O[2] = (Z[A] + Z[B] + Z[C] + Z[D]) * 0.25;

    bool hit = false;

    // First triangle
    hit = Point_in_tri(p->x, p->y, p->z,
                       O[0], O[1], O[2],
                       X[A], Y[A], Z[A],
                       X[B], Y[B], Z[B]);
    if (hit) return true;

    // Second triangle
    hit = Point_in_tri(p->x, p->y, p->z,
                       O[0], O[1], O[2],
                       X[B], Y[B], Z[B],
                       X[C], Y[C], Z[C]);
    if (hit) return true;


    // Third triangle
    hit = Point_in_tri(p->x, p->y, p->z,
                       O[0], O[1], O[2],
                       X[C], Y[C], Z[C],
                       X[D], Y[D], Z[D]);
    if (hit) return true;

    // Fourth triangle
    hit = Point_in_tri(p->x, p->y, p->z,
                       O[0], O[1], O[2],
                       X[D], Y[D], Z[D],
                       X[A], Y[A], Z[A]);
    if (hit) return true;

    return false;
}

bool DynMesh::point_in_face(const Point *p, E_Int fid) const
{
    if (F[fid].size() == 4) return point_in_quad(p, fid);
    return point_in_tri(p, fid);
}

void DynMesh::smooth_skin_ref_data(SkinGraph *skin_graph, E_Int *fdat)
{
    assert(0);
}
