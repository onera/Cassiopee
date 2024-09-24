#include "Mesh.h"
#include "common/mem.h"

// TODO(Imad): sorting routines

struct Point {
    E_Float x, y, z;
};

void points_write(const char *fname, const std::vector<Point> &P)
{
    FILE *fh = fopen(fname, "w");
    assert(fh);
    fprintf(fh, "POINTS\n");
    fprintf(fh, "%zu\n", P.size());
    for (auto p : P) fprintf(fh, "%f %f %f\n", p.x, p.y, p.z);
    fclose(fh);
}

void point_write(E_Float px, E_Float py, E_Float pz)
{
    FILE *fh = fopen("point", "w");
    assert(fh);
    fprintf(fh, "POINTS\n");
    fprintf(fh, "1\n");
    fprintf(fh, "%f %f %f\n", px, py, pz);
    fclose(fh);
}

void point_write(const Point p)
{
    return point_write(p.x, p.y, p.z);
}

struct ArrayI {
    E_Int count;
    E_Int *ptr;
};

void Mesh_extract_points_from_ftag(Mesh *M, ArrayI *pids)
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

struct FaceSort {
    E_Int fid;
    E_Float fc[3];
    E_Float UX, UY, UZ;
    E_Float VX, VY, VZ;
    E_Float UU, VV, UV;
    E_Float inv_denom;
    E_Float xa, ya, za;
};

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

void FaceSort_compute_data(const Mesh *M, FaceSort *mfaces, E_Int mcount)
{
    for (E_Int i = 0; i < mcount; i++) {
        FaceSort *face = &mfaces[i];

        E_Int tid = face->fid;
        assert(M->ftype[tid] == TRI);

        E_Int *tri = Mesh_get_face(M, tid);
        E_Int A = tri[0], B = tri[2], C = tri[4];

        face->UX = (M->X[B] - M->X[A]);
        face->UY = (M->Y[B] - M->Y[A]);
        face->UZ = (M->Z[B] - M->Z[A]);

        face->VX = (M->X[C] - M->X[A]);
        face->VY = (M->Y[C] - M->Y[A]);
        face->VZ = (M->Z[C] - M->Z[A]);

        face->UU = face->UX*face->UX + face->UY*face->UY + face->UZ*face->UZ;
        face->VV = face->VX*face->VX + face->VY*face->VY + face->VZ*face->VZ;
        face->UV = face->UX*face->VX + face->UY*face->VY + face->UZ*face->VZ;

        face->inv_denom = face->UU*face->VV - face->UV*face->UV;

        assert(face->inv_denom != 0.0);

        face->inv_denom = 1.0 / face->inv_denom;

        face->fc[0] = (M->X[A] + M->X[B] + M->X[C]) / 3.0;
        face->fc[1] = (M->Y[A] + M->Y[B] + M->Y[C]) / 3.0;
        face->fc[2] = (M->Z[A] + M->Z[B] + M->Z[C]) / 3.0;

        // Store A
        face->xa = M->X[A];
        face->ya = M->Y[A];
        face->za = M->Z[A];
    }
}

struct Box3 {
    E_Float xmin, ymin, zmin;
    E_Float xmax, ymax, zmax;
};

inline
bool Box3_in_Box3(const Box3 small, const Box3 big)
{
    return (big.xmin <= small.xmin) && (big.xmax >= small.xmax) &&
           (big.ymin <= small.ymin) && (big.ymax >= small.ymax) &&
           (big.zmin <= small.zmin) && (big.zmax >= small.zmax);
}

struct BVH_node {
    Box3 box;
    E_Int start, end;
    BVH_node *left;
    BVH_node *right;
};

Box3 Box3_make
(
    const Mesh *M,
    const FaceSort *mfaces,
    E_Int start, E_Int end
)
{
    E_Float xmin, ymin, zmin, xmax, ymax, zmax;
    xmin = ymin = zmin = FLT_MAX;
    xmax = ymax = zmax = -FLT_MAX;

    for (E_Int i = start; i < end; i++) {
        E_Int fid = mfaces[i].fid;
        E_Int *face = Mesh_get_face(M, fid);
        E_Int *frange = Mesh_get_frange(M, fid);
        assert(M->fstride[fid] == 3);
        assert(M->ftype[fid] == TRI);
        for (E_Int j = 0; j < M->fstride[fid]; j++) {
            E_Int *pn = face + 2*j;
            assert(frange[j] == 1);
            for (E_Int k = 0; k < frange[j]; k++) {
                E_Int pid = pn[k];
                if (M->X[pid] < xmin) xmin = M->X[pid];
                if (M->Y[pid] < ymin) ymin = M->Y[pid];
                if (M->Z[pid] < zmin) zmin = M->Z[pid];

                if (M->X[pid] > xmax) xmax = M->X[pid];
                if (M->Y[pid] > ymax) ymax = M->Y[pid];
                if (M->Z[pid] > zmax) zmax = M->Z[pid];
            }
        }
    }

    // Safety
    E_Float dx = (xmax - xmin) * 0.01;
    E_Float dy = (ymax - ymin) * 0.01;
    E_Float dz = (zmax - zmin) * 0.01;
    xmin -= dx;
    ymin -= dy;
    zmin -= dz;
    xmax += dx;
    ymax += dy;
    zmax += dz;

    return {xmin, ymin, zmin, xmax, ymax, zmax};
}

//static int idx = 0;

BVH_node *BVH_create_node(const Box3 *box, E_Int start, E_Int end,
    BVH_node *left, BVH_node *right)
{
    BVH_node *node = (BVH_node *)XMALLOC(sizeof(BVH_node));
    node->box = *box;
    node->start = start;
    node->end = end;
    node->left = left;
    node->right = right;
    //printf("%d -> %f %f %f %f %f %f\n", idx, box->xmin, box->ymin, box->zmin,
    //    box->xmax, box->ymax, box->zmax);
    //idx++;
    return node;
}

#define MAX_FACES_PER_BVH_LEAF 10

void Box3_clamp(const Box3 *parent, Box3 *child)
{
    child->xmin = std::max(parent->xmin, child->xmin);
    child->ymin = std::max(parent->ymin, child->ymin);
    child->zmin = std::max(parent->zmin, child->zmin);
    child->xmax = std::min(parent->xmax, child->xmax);
    child->ymax = std::min(parent->ymax, child->ymax);
    child->zmax = std::min(parent->zmax, child->zmax);
}

BVH_node *BVH_make(const Mesh *M, FaceSort *mfaces, E_Int start, E_Int end,
    const Box3 *parent_box)
{
    Box3 box = Box3_make(M, mfaces, start, end);
    Box3_clamp(parent_box, &box);
    assert(Box3_in_Box3(box, *parent_box));

    E_Int count = end - start;
    if (count <= MAX_FACES_PER_BVH_LEAF) {
        return BVH_create_node(&box, start, end, NULL, NULL);
    }

    E_Float dx = box.xmax - box.xmin;
    E_Float dy = box.ymax - box.ymin;
    E_Float dz = box.zmax - box.zmin;

    E_Int dim = -1;

    if (dx >= dy && dx >= dz) {
        dim = 0;
    } else if (dy >= dz) {
        dim = 1;
    } else {
        dim = 2;
    }

    std::sort(mfaces + start, mfaces + end,
        [&](const FaceSort &fi, const FaceSort &fj)
        {
            return fi.fc[dim] < fj.fc[dim];
        });
    
    E_Int mid = start + count/2;

    BVH_node *left = BVH_make(M, mfaces, start, mid, &box);
    BVH_node *right = BVH_make(M, mfaces, mid, end, &box);

    assert(Box3_in_Box3(left->box, box));
    assert(Box3_in_Box3(right->box, box));

    return BVH_create_node(&box, start, end, left, right);
}

#define MAX_FACES_PER_POINT 8

struct PointFaces {
    E_Int count;
    E_Int ptr[MAX_FACES_PER_POINT];
};

#define TOL 1e-12

bool Point_in_FaceSort(const Point *p, const FaceSort *f)
{
    E_Float DX = p->x - f->xa;
    E_Float DY = p->y - f->ya;
    E_Float DZ = p->z - f->za;
    
    E_Float d20 = DX*f->UX + DY*f->UY + DZ*f->UZ;
    E_Float d21 = DX*f->VX + DY*f->VY + DZ*f->VZ;

    E_Float u = (f->VV*d20 - f->UV*d21) * f->inv_denom;
    if (u < -TOL || u > 1.0 + TOL) return false;

    E_Float v = (f->UU*d21 - f->UV*d20) * f->inv_denom;
    if (v < -TOL || v > 1.0 + TOL) return false;

    E_Float w = 1.0 - (u + v);
    if (w < -TOL || w > 1.0 + TOL) return false;

    return true;
}

bool Point_in_Box3D(const Point *p, const Box3 *box)
{
    return (p->x >= box->xmin) && (p->x <= box->xmax) &&
           (p->y >= box->ymin) && (p->y <= box->ymax) &&
           (p->z >= box->zmin) && (p->z <= box->zmax);
}

#define MAX_POINTS_PER_FACE 1

void extract_faces_by_threshold
(
    const PointFaces *sploc, E_Int spcount,
    const E_Int *skin, E_Int mcount,
    const E_Int threshold,
    ArrayI *faces
)
{
    E_Int *ftag = (E_Int *)XMALLOC(mcount * sizeof(E_Int));
    memset(ftag, 0, mcount * sizeof(E_Int));

    for (E_Int i = 0; i < spcount; i++) {
        const PointFaces *pfaces = &sploc[i];
        for (E_Int j = 0; j < pfaces->count; j++)
            ftag[pfaces->ptr[j]]++;
    }

    faces->count = 0;

    for (E_Int i = 0; i < mcount; i++) {
        if (ftag[i] > threshold)
            faces->count++;
    }

    faces->ptr = (E_Int *)XMALLOC(faces->count * sizeof(E_Int));
    E_Int *ptr = faces->ptr;

    for (E_Int i = 0; i < mcount; i++) {
        if (ftag[i] > threshold)
            *ptr++ = i;//skin[i];
    }
}

struct SkinGraph {
    E_Int nf;
    E_Int *skin;
    E_Int *xadj;
    E_Int *fpts;
    E_Int *fnei;
};

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

struct EdgeNode {
    E_Int p, q;
    E_Int i, j;
    E_Int posi, posj;
    EdgeNode *next;
};

EdgeNode *make_edge_node(E_Int p, E_Int q, E_Int i, E_Int posi)
{
    EdgeNode *node = (EdgeNode *)XMALLOC(sizeof(EdgeNode));
    node->p = p < q ? p : q;
    node->q = p < q ? q : p;
    node->i = i;
    node->posi = posi;
    node->j = -1;
    node->posj = -1;
    node->next = NULL;
    return node;
}

EdgeNode *find_edge_node(EdgeNode **ht, E_Int hsize, E_Int p, E_Int q)
{
    E_Int p_ = p < q ? p : q;
    E_Int q_ = p < q ? q : p;
    E_Int bucket = p_ % hsize;
    EdgeNode *current = ht[bucket];

    while (current) {
        E_Int P = current->p, Q = current->q;
        if (P == p_ && Q == q_) {
            return current;
        }
        current = current->next;
    }

    return NULL;
}

void insert_edge_node(EdgeNode *node, const E_Int hsize, EdgeNode **ht)
{
    assert(node->p < node->q);
    E_Int bucket = node->p % hsize;
    EdgeNode *current = ht[bucket];

    if (current) {
        EdgeNode *tmp = current;
        ht[bucket] = node;
        node->next = tmp;
    } else {
        ht[bucket] = node;
    }
}

void Mesh_make_skin_neighbours(const Mesh *M, SkinGraph *skin_graph)
{
    E_Int nf = skin_graph->nf;
    const E_Int *xadj = skin_graph->xadj;
    const E_Int *fpts = skin_graph->fpts;

    skin_graph->fnei = (E_Int *)XMALLOC(xadj[nf] * sizeof(E_Int));
    E_Int *fnei = skin_graph->fnei;
    memset(fnei, -1, xadj[nf] * sizeof(E_Int));

    EdgeNode **ht = (EdgeNode **)XMALLOC(nf * sizeof(EdgeNode *));
    memset(ht, 0, nf * sizeof(EdgeNode *));

    for (E_Int i = 0; i < nf; i++) {
        E_Int start = xadj[i];
        E_Int np = xadj[i+1] - start;
        const E_Int *pn = &fpts[start];
        for (E_Int j = 0; j < np; j++) {
            E_Int p = pn[j];
            E_Int q = pn[(j+1)%np];
            EdgeNode *node = find_edge_node(ht, nf, p, q);
            if (node) {
                assert(node->i    != -1);
                assert(node->posi != -1);
                assert(node->j    == -1);
                assert(node->posj == -1);
                node->j = i; 
                node->posj = j;
            } else {
                node = make_edge_node(p, q, i, j);
                insert_edge_node(node, nf, ht);
            }
        } 
    }

    for (E_Int i = 0; i < nf; i++) {
        EdgeNode *node = ht[i];
        while (node) {
            E_Int pi = xadj[node->i] + node->posi;
            assert(fnei[pi] == -1);
            fnei[pi] = node->j;
            
            E_Int pj = xadj[node->j] + node->posj;
            assert(fnei[pj] == -1);
            fnei[pj] = node->i;

            node = node->next;
        }
    }
}

void Mesh_make_skin_graph(const Mesh *M, SkinGraph *skin_graph)
{
    Mesh_extract_skin(M, &skin_graph->nf, &skin_graph->skin);
    Mesh_make_skin_connectivity(M, skin_graph);
    Mesh_make_skin_neighbours(M, skin_graph);
}

struct Vec3f {
    E_Float x, y, z;
};

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

Box3 Box3_make
(
    const Mesh *M,
    const E_Int *skin,
    const E_Int *indices,
    E_Int start, E_Int end
)
{
    E_Float xmin, ymin, zmin, xmax, ymax, zmax;
    xmin = ymin = zmin = FLT_MAX;
    xmax = ymax = zmax = -FLT_MAX;

    for (E_Int i = start; i < end; i++) {
        E_Int fid = skin[indices[i]];
        E_Int *face = Mesh_get_face(M, fid);
        E_Int *frange = Mesh_get_frange(M, fid);

        for (E_Int j = 0; j < M->fstride[fid]; j++) {
            E_Int *pn = face + 2*j;
            for (E_Int k = 0; k < frange[j]; k++) {
                E_Int pid = pn[k];
                if (M->X[pid] < xmin) xmin = M->X[pid];
                if (M->Y[pid] < ymin) ymin = M->Y[pid];
                if (M->Z[pid] < zmin) zmin = M->Z[pid];

                if (M->X[pid] > xmax) xmax = M->X[pid];
                if (M->Y[pid] > ymax) ymax = M->Y[pid];
                if (M->Z[pid] > zmax) zmax = M->Z[pid];
            }
        }
    }

    // Safety
    E_Float dx = (xmax - xmin) * 0.01;
    E_Float dy = (ymax - ymin) * 0.01;
    E_Float dz = (zmax - zmin) * 0.01;
    xmin -= dx;
    ymin -= dy;
    zmin -= dz;
    xmax += dx;
    ymax += dy;
    zmax += dz;

    return {xmin, ymin, zmin, xmax, ymax, zmax};
}

BVH_node *BVH_make(const Mesh *M, const E_Int *skin, const Vec3f *fc,
    E_Int *indices, E_Int start, E_Int end, const Box3 *parent_box)
{
    Box3 box = Box3_make(M, skin, indices, start, end);
    Box3_clamp(parent_box, &box);
    assert(Box3_in_Box3(box, *parent_box));

    E_Int count = end - start;
    if (count <= MAX_FACES_PER_BVH_LEAF) {
        return BVH_create_node(&box, start, end, NULL, NULL);
    }

    E_Float dx = box.xmax - box.xmin;
    E_Float dy = box.ymax - box.ymin;
    E_Float dz = box.zmax - box.zmin;

    E_Int dim = -1;

    if (dx >= dy && dx >= dz) {
        dim = 0;
    } else if (dy >= dz) {
        dim = 1;
    } else {
        dim = 2;
    }

    std::sort(indices + start, indices + end,
        [&](const E_Int i, const E_Int j)
        {
            E_Float *fci = (E_Float *)(&fc[i]);
            E_Float *fcj = (E_Float *)(&fc[j]);
            return fci[dim] < fcj[dim];
        });
    
    E_Int mid = start + count/2;

    BVH_node *left  = BVH_make(M, skin, fc, indices, start, mid, &box);
    BVH_node *right = BVH_make(M, skin, fc, indices, mid, end, &box);

    assert(Box3_in_Box3(left->box, box));
    assert(Box3_in_Box3(right->box, box));

    return BVH_create_node(&box, start, end, left, right);
}

bool point_in_tri(E_Float px, E_Float py, E_Float pz,
    E_Float ax, E_Float ay, E_Float az,
    E_Float bx, E_Float by, E_Float bz,
    E_Float cx, E_Float cy, E_Float cz)
{
    // Normal vector to the plane
    E_Float Y[3] = {bx-ax, by-ay, bz-az};
    E_Float Z[3] = {cx-ax, cy-ay, cz-az};
    E_Float N[3];
    K_MATH::cross(Y, Z, N);

    E_Float X[3] = {px-ax, py-ay, pz-az};

    E_Float dp = K_MATH::dot(N, X, 3);
    
    // Is the point on the plane?
    if (dp < -TOL || dp > TOL) return 0;

    E_Float x1 = K_MATH::dot(X, Y, 3);
    E_Float y1 = K_MATH::dot(Y, Y, 3);
    E_Float z1 = K_MATH::dot(Z, Y, 3);
    E_Float x2 = K_MATH::dot(X, Z, 3);
    E_Float y2 = K_MATH::dot(Y, Z, 3);
    E_Float z2 = K_MATH::dot(Z, Z, 3);

    E_Float u = (x1*z2 - x2*z1) / (y1*z2 - y2*z1);
    if (u < -TOL || u > 1 + TOL) return false;

    E_Float v = (-x1*y2 + x2*y1) / (y1*z2 - y2*z1);
    if (v < -TOL || v > 1 + TOL) return false;

    E_Float w = 1 - u - v;
    if (w < -TOL || w > 1 + TOL) return false;

    return true;
}

bool Mesh_point_in_tri(const Mesh *M, const Point *p, E_Int tid)
{
    const E_Int *face = Mesh_get_face(M, tid);
    E_Int A = face[0], B = face[2], C = face[4];
    return point_in_tri(p->x, p->y, p->z,
                        M->X[A], M->Y[A], M->Z[A],
                        M->X[B], M->Y[B], M->Z[B],
                        M->X[C], M->Y[C], M->Z[C]);
}

bool Mesh_point_in_quad(const Mesh *M, const Point *p, E_Int qid)
{
    // TODO(Imad): maybe compute face centers once in pre-pass
    // Star the quad into 4 triangles
    E_Float O[3] = {0.0, 0.0, 0.0};
    const E_Int *face = Mesh_get_face(M, qid);
    E_Int A = face[0], B = face[2], C = face[4], D = face[6];
    O[0] = (M->X[A] + M->X[B] + M->X[C] + M->X[D]) * 0.25;
    O[1] = (M->Y[A] + M->Y[B] + M->Y[C] + M->Y[D]) * 0.25;
    O[2] = (M->Z[A] + M->Z[B] + M->Z[C] + M->Z[D]) * 0.25;

    bool hit = false;

    // First triangle
    hit = point_in_tri(p->x, p->y, p->z,
                       O[0], O[1], O[2],
                       M->X[A], M->Y[A], M->Z[A],
                       M->X[B], M->Y[B], M->Z[B]);
    if (hit) return true;

    // Second triangle
    hit = point_in_tri(p->x, p->y, p->z,
                       O[0], O[1], O[2],
                       M->X[B], M->Y[B], M->Z[B],
                       M->X[C], M->Y[C], M->Z[C]);
    if (hit) return true;


    // Third triangle
    hit = point_in_tri(p->x, p->y, p->z,
                       O[0], O[1], O[2],
                       M->X[C], M->Y[C], M->Z[C],
                       M->X[D], M->Y[D], M->Z[D]);
    if (hit) return true;

    // Fourth triangle
    hit = point_in_tri(p->x, p->y, p->z,
                       O[0], O[1], O[2],
                       M->X[D], M->Y[D], M->Z[D],
                       M->X[A], M->Y[A], M->Z[A]);
    if (hit) return true;

    return false;
}

bool Mesh_point_in_face(const Mesh *M, const Point *p, E_Int fid)
{
    if (M->ftype[fid] == QUAD) return Mesh_point_in_quad(M, p, fid);
    assert(M->ftype[fid] == TRI);
    return Mesh_point_in_tri(M, p, fid);
}

void BVH_locate_point
(
    const BVH_node *node,
    const Mesh *M,
    const E_Int *skin,
    const E_Int *indices,
    const Point *p,
    PointFaces *pfaces
)
{
    if (node->left == NULL && node->right == NULL) {
        for (E_Int i = node->start; i < node->end; i++) {
            E_Int fid = skin[indices[i]];

            if (Mesh_point_in_face(M, p, fid)) {
                if (pfaces->count >= MAX_FACES_PER_POINT) {
                    fprintf(stderr, 
                        "bvh_locate: MAX_FACES_PER_POINT exceeded!\n");
                    abort();
                }
                pfaces->ptr[pfaces->count++] = i;
            }
        }
        return;
    }

    assert(node->left && node->right);
    assert(Box3_in_Box3(node->left->box, node->box));
    assert(Box3_in_Box3(node->right->box, node->box));

    bool in_box = Point_in_Box3D(p, &node->box);
    
    if (!in_box)
        return;

    BVH_locate_point(node->left, M, skin, indices, p, pfaces);
    BVH_locate_point(node->right, M, skin, indices, p, pfaces);
}

void locate_spoints_in_mskin
(
    const Mesh *S,
    const ArrayI *spoints,
    const Mesh *M,
    const E_Int *skin,
    const E_Int *indices,
    const BVH_node *bvh,
    PointFaces *sploc
)
{
    for (E_Int i = 0; i < spoints->count; i++) {
        E_Int spid = spoints->ptr[i];
        const Point p = {S->X[spid], S->Y[spid], S->Z[spid]};
        PointFaces *pfaces = &sploc[i];
        assert(Point_in_Box3D(&p, &bvh->box));
        BVH_locate_point(bvh, M, skin, indices, &p, pfaces);
        if (pfaces->count == 0) {
            fprintf(stderr, "bvh_locate: failed at point index %d (%d)!\n",
                i, spid);
            point_write(p);
            abort();
        }
    }
}

void smooth_skin_ref_data(const SkinGraph *skin_graph, E_Int *fdat)
{
    E_Int nf = skin_graph->nf;
    E_Int stack_size = 3*nf;
    E_Int *fstack = (E_Int *)XMALLOC(stack_size * sizeof(E_Int));
    memset(fstack, -1, stack_size * sizeof(E_Int));

    E_Int l = 0;

    for (E_Int i = 0; i < nf; i++) {
        if (fdat[i] > 0)
            fstack[l++] = i;
    }

    const E_Int *xadj = skin_graph->xadj;
    const E_Int *fnei = skin_graph->fnei;

    while (--l >= 0) {
        E_Int fid = fstack[l];

        E_Int start = xadj[fid];
        E_Int nneis = xadj[fid+1] - start;
        const E_Int *neis = &fnei[start];
        for (E_Int i = 0; i < nneis; i++) {
            E_Int nei = neis[i];
            E_Int incr_nei = fdat[nei];
            E_Int incr_fid = fdat[fid];
            E_Int diff = abs(incr_nei - incr_fid);
            if (diff <= 1) continue;
            E_Int idx_to_modify = incr_fid > incr_nei ? nei : fid;
            fdat[idx_to_modify] += diff-1;
            l++;
            assert(l < stack_size);
            fstack[l] = idx_to_modify;
        }
    }

    XFREE(fstack);
}

void init_skin_refinement_cells(const SkinGraph *skin_graph,
    const E_Int *fdat, Mesh *M, E_Int *rcount)
{
    E_Int nf = skin_graph->nf;
    E_Int count = 0;

    for (E_Int i = 0; i < nf; i++) {
        if (fdat[i] > 0) {
            E_Int fid = skin_graph->skin[i];
            E_Int own = M->owner[fid];
            // If cref[own] != 0, H18 is impossible, refinement has to be H27.
            assert(M->cref[own] == 0);
            M->cref[own] = 1;
            count++;
        }
    }

    *rcount = count;
}

/*
void Mesh_extract_refinement_cells(const Mesh *M, ArrayI *rcells)
{
    rcells->count = 0;

    for (E_Int cid = 0; cid < M->nc; cid++) {
        rcells->count += (M->cref[cid] > 0);
    }

    rcells->ptr = (E_Int *)XMALLOC(rcells->count * sizeof(E_Int));
    E_Int *ptr = rcells->ptr;

    for (E_Int cid = 0; cid < M->nc; cid++) {
        if (M->cref[cid] > 0)
            *ptr++ = cid;
    }
}
*/

void smooth_cell_refinement_data(Mesh *M)
{
    E_Int nc = M->nc;
    E_Int stack_size = 3*nc;
    E_Int *cstack = (E_Int *)XMALLOC(stack_size * sizeof(E_Int));
    memset(cstack, -1, stack_size * sizeof(E_Int));

    E_Int l = 0;

    for (E_Int cid = 0; cid < nc; cid++) {
        if (M->cref[cid] > 0)
            cstack[l++] = cid;
    }

    while (--l >= 0) {
        E_Int cid = cstack[l];

        E_Int nn, neis[24];
        Mesh_get_cneis(M, cid, nn, neis);

        for (E_Int i = 0; i < nn; i++) {
            E_Int nei = neis[i];
            E_Int incr_nei = M->cref[nei] + M->clevel[nei];
            E_Int incr_cid = M->cref[cid] + M->clevel[cid];
            E_Int diff = abs(incr_nei - incr_cid);
            if (diff <= 1) continue;
            E_Int idx_to_modify = incr_cid > incr_nei ? nei : cid;
            M->cref[idx_to_modify] += 1;
            l++;
            assert(l < stack_size);
            cstack[l] = idx_to_modify;
        }
    }
}

static
void Mesh_set_face_as_cell_bottom(Mesh *M, E_Int fid, E_Int cid)
{
    // Make fid the bottom face
    E_Int *cell = Mesh_get_cell(M, cid);
    E_Int size = 4*M->cstride[cid];
    E_Int pos = Get_pos(fid, cell, size);
    assert(pos != -1);
    assert(pos % 4 == 0);
    Right_shift(cell, pos, size);
    assert(cell[0] == fid);
    E_Int *crange = Mesh_get_crange(M, cid);
    Right_shift(crange, pos/4, M->cstride[cid]);
}

void reorder_cells_for_H18(const SkinGraph *skin_graph, const E_Int *indices,
    const ArrayI *rfaces, Mesh *M)
{
    for (E_Int i = 0; i < rfaces->count; i++) {
        E_Int idx_in_skin = indices[rfaces->ptr[i]];
        E_Int fid = skin_graph->skin[idx_in_skin];
        E_Int own = M->owner[fid];
        assert(M->cref[own] == 1);

        while (M->cref[own] == 1) {
            Mesh_set_face_as_cell_bottom(M, fid, own);
            // Get the top face
            fid = Mesh_get_cell(M, own)[4];
            // Get the top neighbour
            own = Mesh_get_cnei(M, own, fid);
        }
    }
}

void assign_face_refinement_data(Mesh *M)
{
    for (E_Int cid = 0; cid < M->nc; cid++) {
        if (M->cref[cid] == 0) continue;
        assert(M->cref[cid] == 1);

        E_Int *cell = Mesh_get_cell(M, cid);
        E_Int *crange = Mesh_get_crange(M, cid);
        E_Int cstride = M->cstride[cid];
        E_Int clvl = M->clevel[cid];

        for (E_Int i = 0; i < cstride; i++) {
            E_Int *pf = cell + 4*i;

            for (E_Int j = 0; j < crange[i]; j++) {
                E_Int face = pf[j];
                if (M->fref[face] == 1) continue;

                E_Int flvl = M->flevel[face];

                assert(flvl >= clvl);

                if (flvl == clvl) M->fref[face] = 1;
            }
        }
    }
}

// TODO(Imad): isotropic resizing for now

void Mesh_isolate_refinement_entities(Mesh *M, ArrayI *rcells, ArrayI *rfaces)
{
    rcells->count = 0;
    for (E_Int cid = 0; cid < M->nc; cid++) {
        rcells->count += (M->cref[cid] == 1);
    }
    rcells->ptr = (E_Int *)XMALLOC(rcells->count * sizeof(E_Int));
    E_Int *ptr = rcells->ptr;
    for (E_Int cid = 0; cid < M->nc; cid++) {
        if (M->cref[cid] > 0) {
            assert(M->cref[cid] == 1);
            *ptr++ = cid;
        }
    }

    rfaces->count = 0;
    for (E_Int fid = 0; fid < M->nf; fid++) {
        rfaces->count += (M->fref[fid] == 1);
    }
    rfaces->ptr = (E_Int *)XMALLOC(rfaces->count * sizeof(E_Int));
    ptr = rfaces->ptr;
    for (E_Int fid = 0; fid < M->nf; fid++) {
        if (M->fref[fid] > 0) {
            assert(M->fref[fid] == 1);
            *ptr++ = fid;
        }
    }
}

void Mesh_resize(Mesh *M, const ArrayI *rcells, const ArrayI *rfaces)
{
    // 7 new cells per refined cells
    E_Int cell_incr = rcells->count * 7;
    // 3 new faces per refined face + 12 new faces per refined cell
    E_Int face_incr = rfaces->count * 3 + rcells->count * 12;
    // Estimate point increase
    E_Int point_incr = 2 * face_incr;

    E_Int new_nc = M->nc + cell_incr;
    E_Int new_nf = M->nf + face_incr;
    E_Int new_np = M->np + point_incr;

    M->X = (E_Float *)XRESIZE(M->X, new_np * sizeof(E_Float));
    M->Y = (E_Float *)XRESIZE(M->Y, new_np * sizeof(E_Float));
    M->Z = (E_Float *)XRESIZE(M->Z, new_np * sizeof(E_Float));

    Mesh_resize_face_data(M, new_nf);
    Mesh_resize_cell_data(M, new_nc);
}

#define DIR_ISO 0
#define DIR_X 1
#define DIR_Y 2

void prepare_cell_ordering_and_face_refinement_patterns(Mesh *M,
    const ArrayI *rcells, const ArrayI *rfaces)
{
    for (E_Int i = 0; i < rcells->count; i++) {
        E_Int cid = rcells->ptr[i];
        assert(M->cref[cid] == 1);
        if (M->clevel[cid] == 0) {
            H18_reorder(cid, M);
            assert(check_canon_hexa(cid, M) == 0); 
        }
    }

    for (E_Int i = 0; i < rfaces->count; i++) {
        E_Int fid = rfaces->ptr[i];
        if (M->flevel[fid] > 0) continue;

        E_Int cid = M->owner[fid];
        if (M->cref[cid] == 0)
            cid = M->neigh[fid]; 

        assert(cid != -1);
        assert(M->cref[cid] == 1);
        
        E_Int *cell = Mesh_get_cell(M, cid);
        E_Int *crange = Mesh_get_crange(M, cid);

        E_Int pos = Get_pos(fid, cell, 4*M->cstride[cid]);
        E_Int side = pos / 4;
        E_Int *face = Mesh_get_face(M, fid);

        if (side == 0 || side == 1) {
            M->fpattern[fid] = DIR_ISO;
        } else {

            // Reconstruct bottom
            E_Int map[4];
            E_Int bot = cell[0];
            E_Int N0 = Mesh_get_face(M, bot)[0];
            reconstruct_quad(M, cid, cell, crange[0], normalIn_H[0], N0, map); 

            M->fpattern[fid] = DIR_ISO;
            
            E_Int i0;
            E_Int reorient;

            if (side == 2) {

                E_Int np, lpts[8];
                for (E_Int j = 0; j < 8; j++) lpts[j] = -1;
                Mesh_get_fpoints(M, fid, np, lpts);
            
                // Must share map[0] && map[3]
                i0 = Get_pos(map[3], face, 8);
                assert(i0 != -1);
                i0 = Get_pos(map[0], face, 8);
                assert(i0 != -1);
                reorient = Mesh_get_reorient(M, fid, cid, normalIn_H[side]);

            } else if (side == 3) {

                // Must share map[1] && map[2]
                i0 = Get_pos(map[2], face, 8);
                assert(i0 != -1);
                i0 = Get_pos(map[1], face, 8);
                assert(i0 != -1);
                reorient = Mesh_get_reorient(M, fid, cid, normalIn_H[side]);

            } else if (side == 4) {

                // Must share map[1] && map[0]
                i0 = Get_pos(map[0], face, 8);
                assert(i0 != -1);
                i0 = Get_pos(map[1], face, 8);
                assert(i0 != -1);
                reorient = Mesh_get_reorient(M, fid, cid, normalIn_H[side]);

            } else if (side == 5) {
                
                // Must share map[2] && map[3]
                i0 = Get_pos(map[3], face, 8);
                assert(i0 != -1);
                i0 = Get_pos(map[2], face, 8);
                assert(i0 != -1);
                reorient = Mesh_get_reorient(M, fid, cid, normalIn_H[side]);

            } else {
                assert(0);
            }


            i0 /= 2;

            assert(i0 == 0 || i0 == 1 || i0 == 2 || i0 == 3);

            if (reorient == 0) {
                if (i0 == 0 || i0 == 2) M->fpattern[fid] = DIR_X;
                else M->fpattern[fid] = DIR_Y;
            } else {
                if (i0 == 0 || i0 == 2) M->fpattern[fid] = DIR_Y;
                else M->fpattern[fid] = DIR_X;
            }

            assert(M->fpattern[fid] != DIR_ISO);
            
            /*
            bool swap_1 = (reorient == 0) && (i0 == 1 || i0 == 3);
            bool swap_2 = (reorient == 1) && (i0 == 0 || i0 == 2);

            if (swap_1 || swap_2) M->fpattern[fid] = DIR_Y;
            */
        }
    }
}

E_Int refine_quad_X(E_Int quad, Mesh *M)
{
    assert(M->fref[quad] == 1);
    E_Int NODES[8];

    E_Int *fpts = Mesh_get_face(M, quad);
    E_Int *frange = Mesh_get_frange(M, quad);

    // BOT
    NODES[0] = fpts[0];
    if (frange[0] == 2) {
        NODES[4] = fpts[1];
    } else {
        E_Int p = fpts[0];
        E_Int q = fpts[2];
        Mesh_refine_or_get_edge_center(M, p, q, NODES[4]);
    }

    // RGT
    NODES[1] = fpts[2];
    NODES[5] = fpts[3];

    // TOP
    NODES[2] = fpts[4];
    if (frange[2] == 2) {
        NODES[6] = fpts[5];
    } else {
        E_Int p = fpts[4];
        E_Int q = fpts[6];
        Mesh_refine_or_get_edge_center(M, p, q, NODES[6]);
    }

    // LFT
    NODES[3] = fpts[6];
    NODES[7] = fpts[7];


    // Second child
    fpts = Mesh_get_face(M, M->nf);
    fpts[0] = NODES[4];
    fpts[1] = -1;
    fpts[2] = NODES[1];
    fpts[3] = NODES[5];
    fpts[4] = NODES[2];
    fpts[5] = -1;
    fpts[6] = NODES[6];
    fpts[7] = -1;

    // First child replaces quad
    fpts = Mesh_get_face(M, quad);
    fpts[0] = NODES[0];
    fpts[1] = -1;
    fpts[2] = NODES[4];
    fpts[3] = -1;
    fpts[4] = NODES[6];
    fpts[5] = -1;
    fpts[6] = NODES[3];
    fpts[7] = NODES[7];

    // Update ranges and strides
    Mesh_update_face_range_and_stride(M, quad, M->nf, 1);

    // Conformize parent cells

    E_Int own = M->owner[quad];

    assert(M->clevel[own] == M->flevel[quad]);

    if (Mesh_conformize_cell_face(M, own, quad, M->nf, 2) != 0) return 1;

    E_Int nei = M->neigh[quad];

    if (nei != -1) {
        assert(M->clevel[nei] == M->flevel[quad]);

        if (Mesh_conformize_cell_face(M, nei, quad, M->nf, 2) != 0) return 1;
    }

    for (E_Int i = 0; i < 1; i++) {
        M->owner[M->nf+i] = own;
        M->neigh[M->nf+i] = nei;
    }

    // Update adaptation info
    M->flevel[quad]++;

    assert(M->fref[quad] == FACE_REFINED);

    for (E_Int i = 0; i < 1; i++) {
        E_Int fid = M->nf + i;

        M->flevel[fid] = M->flevel[quad];
        M->ftype[fid] = M->ftype[quad];

        M->fref[fid] = FACE_NEW;
    }

    M->fchildren[quad] = {quad, M->nf};

    M->fparent[M->nf] = quad;

    M->ftag[M->nf] = M->ftag[quad];

    // Increment face/edge/point count
    M->nf += 1;

    return 0;
}

E_Int refine_quad_Y(E_Int quad, Mesh *M)
{
    assert(M->fref[quad] == 1);
    E_Int NODES[8];

    E_Int *fpts = Mesh_get_face(M, quad);
    E_Int *frange = Mesh_get_frange(M, quad);

    // BOT
    NODES[0] = fpts[0];
    NODES[4] = fpts[1];

    // RGT
    NODES[1] = fpts[2];
    if (frange[1] == 2) {
        NODES[5] = fpts[3];
    } else {
        E_Int p = fpts[2];
        E_Int q = fpts[4];
        Mesh_refine_or_get_edge_center(M, p, q, NODES[5]);
    }

    // TOP
    NODES[2] = fpts[4];
    NODES[6] = fpts[5];
    
    // LFT
    NODES[3] = fpts[6];
    if (frange[3] == 2) {
        NODES[7] = fpts[7];
    } else {
        E_Int p = fpts[6];
        E_Int q = fpts[0];
        Mesh_refine_or_get_edge_center(M, p, q, NODES[7]);
    }

    // Second child
    fpts = Mesh_get_face(M, M->nf);
    fpts[0] = NODES[7];
    fpts[1] = -1;
    fpts[2] = NODES[5];
    fpts[3] = -1;
    fpts[4] = NODES[2];
    fpts[5] = NODES[6];
    fpts[6] = NODES[3];
    fpts[7] = -1;

    // First child replaces quad
    fpts = Mesh_get_face(M, quad);
    fpts[0] = NODES[0];
    fpts[1] = NODES[4];
    fpts[2] = NODES[1];
    fpts[3] = -1;
    fpts[4] = NODES[5];
    fpts[5] = -1;
    fpts[6] = NODES[7];
    fpts[7] = -1;

    // Update ranges and strides
    Mesh_update_face_range_and_stride(M, quad, M->nf, 1);

    // Conformize parent cells

    E_Int own = M->owner[quad];

    assert(M->clevel[own] == M->flevel[quad]);

    if (Mesh_conformize_cell_face(M, own, quad, M->nf, 2) != 0) return 1;

    E_Int nei = M->neigh[quad];

    if (nei != -1) {
        assert(M->clevel[nei] == M->flevel[quad]);

        if (Mesh_conformize_cell_face(M, nei, quad, M->nf, 2) != 0) return 1;
    }

    for (E_Int i = 0; i < 1; i++) {
        M->owner[M->nf+i] = own;
        M->neigh[M->nf+i] = nei;
    }

    // Update adaptation info
    M->flevel[quad]++;

    assert(M->fref[quad] == FACE_REFINED);

    for (E_Int i = 0; i < 1; i++) {
        E_Int fid = M->nf + i;

        M->flevel[fid] = M->flevel[quad];
        M->ftype[fid] = M->ftype[quad];

        M->fref[fid] = FACE_NEW;
    }

    M->fchildren[quad] = {quad, M->nf};

    M->fparent[M->nf] = quad;

    M->ftag[M->nf] = M->ftag[quad];

    // Increment face/edge/point count
    M->nf += 1;

    return 0;
}

E_Int refine_quad_dir(E_Int fid, Mesh *M)
{
    if (M->fpattern[fid] == DIR_X) return refine_quad_X(fid, M);
    if (M->fpattern[fid] == DIR_Y) return refine_quad_Y(fid, M);
    return Q9_refine(fid, M);
}

void Mesh_refine_dir(Mesh *M, ArrayI *ref_cells, ArrayI *ref_faces)
{
    std::set<E_Int> levelset;

    for (E_Int i = 0; i < ref_cells->count; i++) {
        levelset.insert(M->clevel[ref_cells->ptr[i]]);
    }

    for (E_Int i = 0; i < ref_faces->count; i++) {
        levelset.insert(M->flevel[ref_faces->ptr[i]]);
    }

    std::vector<E_Int> levels;
    for (E_Int level : levelset) levels.push_back(level);
    std::sort(levels.begin(), levels.end());

    //std::reverse(ref_cells->ptr, ref_cells->ptr + ref_cells->count);
    //std::reverse(ref_faces->ptr, ref_faces->ptr + ref_faces->count);

    E_Int cells_left = ref_cells->count-1;
    E_Int faces_left = ref_faces->count-1;

    for (E_Int level : levels) {
        while (faces_left >= 0 && M->flevel[ref_faces->ptr[faces_left]] == level) {
            E_Int face = ref_faces->ptr[faces_left];
            faces_left--;
            refine_quad_dir(face, M);
        }

        while (cells_left >= 0 && M->clevel[ref_cells->ptr[cells_left]] == level) {
            E_Int cell = ref_cells->ptr[cells_left];
            cells_left--;
            refine_cell_dir(cell, M);
        }
    }
}

PyObject *K_XCORE::AdaptMesh_AdaptGeom(PyObject *self, PyObject *args)
{
    PyObject *AMESH, *SMESH;

    if (!PYPARSETUPLE_(args, OO_, &AMESH, &SMESH)) {
        RAISE("Wrong input.");
        return NULL;
    }

    if (!PyCapsule_IsValid(AMESH, "AdaptMesh")) {
        RAISE("Bad first AdaptMesh hook.");
        return NULL;
    }

    if (!PyCapsule_IsValid(SMESH, "AdaptMesh")) {
        RAISE("Bad second AdaptMesh hook.");
        return NULL;
    }

    Mesh *M = (Mesh *)PyCapsule_GetPointer(AMESH, "AdaptMesh");
    Mesh *S = (Mesh *)PyCapsule_GetPointer(SMESH, "AdaptMesh");

    if (M->npc > 1) {
        RAISE("AdaptGeom is sequential.");
        return NULL;
    }

    // Refine M volumetric wrt to S tagged faces point cloud
    // Refine S surfacic wrt to M tagged faces point cloud

    //E_Int ref_M = 0;
    E_Int ref_S = 0;
    E_Int iter = 0;

    do {
        iter++;

        // Extract spoints from tagged sfaces
        ArrayI spoints;
        Mesh_extract_points_from_ftag(S, &spoints);
        
        // We need the skin connectivity graph
        SkinGraph skin_graph;
        Mesh_make_skin_graph(M, &skin_graph);
        printf("Skin: %d faces\n", skin_graph.nf);

        // BVH the skin
        Vec3f *skin_fc = (Vec3f *)XMALLOC(skin_graph.nf * sizeof(Vec3f));
        Mesh_make_face_centers(M, skin_graph.nf, skin_graph.skin, skin_fc);
        E_Int *indices = (E_Int *)XMALLOC(skin_graph.nf * sizeof(E_Int));
        for (E_Int i = 0; i < skin_graph.nf; i++) indices[i] = i;
        const Box3 huge = {-FLT_MAX, -FLT_MAX, -FLT_MAX,
                            FLT_MAX,  FLT_MAX,  FLT_MAX};
        BVH_node *bvh = BVH_make(M, skin_graph.skin, skin_fc, indices, 0,
            skin_graph.nf, &huge);
        puts("BVH constructed");

        // Locate spoints in skin
        PointFaces *sploc =
            (PointFaces *)XMALLOC(spoints.count * sizeof(PointFaces));
        memset(sploc, 0, spoints.count * sizeof(PointFaces));
        locate_spoints_in_mskin
        (
            S, &spoints,
            M, skin_graph.skin, indices,
            bvh,
            sploc
        );
        puts("Points located");

        // Isolate faces that contain more than MAX_POINTS_PER_FACE spoints
        ArrayI rfaces;
        extract_faces_by_threshold
        (
            sploc, spoints.count,
            skin_graph.skin, skin_graph.nf,
            MAX_POINTS_PER_FACE,
            &rfaces
        );
        puts("Refinement faces isolated");
        printf("Refinement faces: %d\n", rfaces.count);

        // Smooth face refinement data
        E_Int *fdat = (E_Int *)XMALLOC(skin_graph.nf * sizeof(E_Int));
        memset(fdat, 0, skin_graph.nf * sizeof(E_Int));
        for (E_Int i = 0; i < rfaces.count; i++) {
            E_Int idx_in_skin = indices[rfaces.ptr[i]];
            E_Int fid = skin_graph.skin[idx_in_skin];
            fdat[idx_in_skin] = M->flevel[fid] + 1;
        }
        smooth_skin_ref_data(&skin_graph, fdat);
        E_Int smooth_nfref = 0;
        for (E_Int i = 0; i < rfaces.count; i++) {
            E_Int idx_in_skin = indices[rfaces.ptr[i]];
            E_Int fid = skin_graph.skin[idx_in_skin];
            fdat[idx_in_skin] -= M->flevel[fid];
            if (fdat[idx_in_skin] > 0) smooth_nfref++;
        }
        puts("Face refinement smoothed out");
        printf("Smooth refinement face count: %d\n", smooth_nfref);

        /*
        npy_intp dims[2];
        dims[1] = 1;
        dims[0] = (npy_intp)smooth_nfref;
        PyArrayObject *FACES = (PyArrayObject *)PyArray_SimpleNew(1, dims, E_NPY_INT);

        E_Int *pf = (E_Int *)PyArray_DATA(FACES);
        E_Int *ptr = pf;
        for (E_Int i = 0; i < skin_graph.nf; i++) {
            if (fdat[i] > 0) {
                *ptr++ = skin_graph.skin[i]+1;
            }
        }
        return (PyObject *)FACES;
        */

        // Allocate
        M->cref = (E_Int *)XRESIZE(M->cref, M->nc * sizeof(E_Int));
        memset(M->cref, 0, M->nc * sizeof(E_Int));

        // Cells
        E_Int ref_cell_count;
        init_skin_refinement_cells(&skin_graph, fdat, M, &ref_cell_count);
        printf("Refinement cells: %d\n", ref_cell_count);

        // Smooth cell refinement data
        smooth_cell_refinement_data(M);
        puts("Cell refinement smoothed out");

        // Assign refinement data
        M->fref = (E_Int *)XRESIZE(M->fref, M->nf * sizeof(E_Int));
        memset(M->fref, 0, M->nf * sizeof(E_Int));
        assign_face_refinement_data(M);

        // Isolate cells/faces to be refined
        ArrayI ref_cells, ref_faces;
        Mesh_isolate_refinement_entities(M, &ref_cells, &ref_faces);
        printf("Refinement cells: %d\n", ref_cells.count);
        printf("Refinement faces: %d\n", ref_faces.count);

        // Setup bottom/top face chain
        reorder_cells_for_H18(&skin_graph, indices, &rfaces, M);

        // Prepare cell ordering and face refinement direction
        M->fpattern = (E_Int *)XMALLOC(M->nf * sizeof(E_Int));
        memset(M->fpattern, -1, M->nf * sizeof(E_Int));
        prepare_cell_ordering_and_face_refinement_patterns(M, &ref_cells,
            &ref_faces);
        
        // Resize for refinement
        Mesh_resize(M, &ref_cells, &ref_faces);
        puts("Mesh resized for refinement");

        // Sort entities by refinement level
        std::sort(ref_cells.ptr, ref_cells.ptr + ref_cells.count,
            [&] (E_Int i, E_Int j) { return M->clevel[i] > M->clevel[j]; });
        puts("Refinement cells sorted");
        std::sort(ref_faces.ptr, ref_faces.ptr + ref_faces.count,
            [&] (E_Int i, E_Int j) { return M->flevel[i] > M->flevel[j]; });
        puts("Refinement faces sorted");
        
        // Refine
        printf("Cells before refinement: %d\n", M->nc);
        printf("Faces before refinement: %d\n", M->nf);
        Mesh_refine_dir(M, &ref_cells, &ref_faces);

        printf("Cells after refinement: %d\n", M->nc);
        printf("Faces after refinement: %d\n", M->nf);


        

        // FREE

    } while (ref_S);

    
    Mesh_conformize_face_edge(M);    

    return Py_None;
}
