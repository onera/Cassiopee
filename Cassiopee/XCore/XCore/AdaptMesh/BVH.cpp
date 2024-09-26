#include "BVH.h"
#include "Mesh.h"
#include "FaceSort.h"
#include "Point.h"
#include "common/mem.h"
#include "Vec.h"
#include "DynMesh.h"

#define MAX_FACES_PER_BVH_LEAF 10

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

void BVH_free(BVH_node *node)
{
    if (node == NULL) return;

    BVH_free(node->left);
    BVH_free(node->right);
    
    XFREE(node);
};

/* DynMesh */

BVH_node *BVH_make(const DynMesh *M, FaceSort *mfaces, E_Int start, E_Int end,
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

BVH_node *BVH_make(const DynMesh *M, const E_Int *skin, const Vec3f *fc,
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

void BVH_locate_point
(
    const BVH_node *node,
    const DynMesh *M,
    const E_Int *skin,
    const E_Int *indices,
    const Point *p,
    PointFaces *pfaces
)
{
    if (node->left == NULL && node->right == NULL) {
        for (E_Int i = node->start; i < node->end; i++) {
            E_Int fid = skin[indices[i]];

            if (M->point_in_face(p, fid)) {
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