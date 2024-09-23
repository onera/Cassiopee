#include "Mesh.h"
#include "common/mem.h"

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

void BVH_locate_point
(
    const BVH_node *node,
    const FaceSort *mfaces,
    const Point *p,
    PointFaces *pfaces
)
{
    if (node->left == NULL && node->right == NULL) {
        for (E_Int i = node->start; i < node->end; i++) {
            const FaceSort *mface = &mfaces[i];

            if (Point_in_FaceSort(p, mface)) {
                if (pfaces->count >= MAX_FACES_PER_POINT) {
                    fprintf(stderr, 
                        "bvh_locate: MAX_FACES_PER_POINT exceeded!\n");
                    abort();
                }
                pfaces->ptr[pfaces->count++] = i; //mface->fid;
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

    BVH_locate_point(node->left, mfaces, p, pfaces);
    BVH_locate_point(node->right, mfaces, p, pfaces);
}

void locate_spoints_in_mtris
(
    const Mesh *S,
    const ArrayI *spoints,
    const FaceSort *mfaces,
    const BVH_node *bvh,
    PointFaces *sploc
)
{
    for (E_Int i = 0; i < spoints->count; i++) {
        E_Int spid = spoints->ptr[i];
        const Point p = {S->X[spid], S->Y[spid], S->Z[spid]};
        PointFaces *pfaces = &sploc[i];
        assert(Point_in_Box3D(&p, &bvh->box));
        BVH_locate_point(bvh, mfaces, &p, pfaces);
        if (pfaces->count == 0) {
            fprintf(stderr, "bvh_locate: failed at point index %d (%d)!\n",
                i, spid);
            point_write(p);
            abort();
        }
    }
}

#define MAX_POINTS_PER_FACE 1

void extract_faces_by_threshold
(
    const PointFaces *sploc, E_Int spcount,
    const FaceSort *mtris, E_Int mcount,
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
            *ptr++ = mtris[i].fid;
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

void Mesh_make_skin_graph(Mesh *M, SkinGraph *skin_graph)
{
    Mesh_extract_skin(M, &skin_graph->nf, &skin_graph->skin);
    Mesh_make_skin_connectivity(M, skin_graph);
    Mesh_make_skin_neighbours(M, skin_graph);
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
        
        /*
        // Spoints must belong to one of the tagged Mfaces
        FaceSort *mfaces = NULL;
        E_Int mcount = 0;
        Mesh_extract_faces_from_ftag(M, &mfaces, &mcount);
        printf("Mfaces: %d\n", mcount);
        assert(mfaces);

        // Compute their centroids
        FaceSort_compute_data(M, mfaces, mcount);

        // Build a BVH for mfaces
        const Box3 huge = {-FLT_MAX, -FLT_MAX, -FLT_MAX,
                            FLT_MAX,  FLT_MAX,  FLT_MAX};
        BVH_node *bvh = BVH_make(M, mfaces, 0, mcount, &huge);

        // Locate spoints in mfaces
        PointFaces *sploc =
            (PointFaces *)XMALLOC(spoints.count * sizeof(PointFaces));
        memset(sploc, 0, spoints.count * sizeof(PointFaces));
        locate_spoints_in_mtris
        (
            S, &spoints,
            mfaces, bvh,
            sploc
        );

        // Isolate faces that contain more than MAX_POINTS_PER_FACE spoints
        ArrayI rfaces;
        extract_faces_by_threshold
        (
            sploc, spoints.count,
            mfaces, mcount,
            MAX_POINTS_PER_FACE,
            &rfaces
        );

        printf("Refinement faces: %d\n", rfaces.count);
        */
        
        // We need the skin connectivity graph
        SkinGraph skin_graph;
        Mesh_make_skin_graph(M, &skin_graph);
        printf("Skin: %d faces\n", skin_graph.nf);




        // FREE

    } while (ref_S);


    return Py_None;
}
