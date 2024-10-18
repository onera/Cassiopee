#include "smesh.h"
#include "AABB.h"
#include "BVH.h"

AABB Smesh::make_AABB(E_Int start, E_Int end)
{
    E_Float xmin, ymin, zmin, xmax, ymax, zmax;
    xmin = ymin = zmin = EFLOATMAX;
    xmax = ymax = zmax = EFLOATMIN;

    for (E_Int i = start; i < end; i++) {
        E_Int fid = bvh_fids[i];
        const auto &pn = F[fid];
        for (E_Int p : pn) {
            if (X[p] < xmin) xmin = X[p];
            if (Y[p] < ymin) ymin = Y[p];
            if (Z[p] < zmin) zmin = Z[p];
            if (X[p] > xmax) xmax = X[p];
            if (Y[p] > ymax) ymax = Y[p];
            if (Z[p] > zmax) zmax = Z[p]; 
        }
    }

    E_Float dx = xmax - xmin;
    E_Float dy = ymax - ymin;
    E_Float dz = zmax - zmin;

    xmin -= dx * 0.01;
    ymin -= dy * 0.01;
    zmin -= dz * 0.01;
    xmax += dx * 0.01;
    ymax += dy * 0.01;
    zmax += dz * 0.01;

    return {xmin, ymin, zmin, xmax, ymax, zmax};
}

BVH_node *Smesh::make_BVH_node(const AABB &box, E_Int start, E_Int end,
    BVH_node *left, BVH_node *right)
{
    BVH_node *node = new BVH_node;
    node->box = box;
    node->start = start;
    node->end = end;
    node->left = left;
    node->right = right;
    return node;
}

BVH_node *Smesh::make_BVH_subtree(E_Int start, E_Int end, const AABB &parent)
{
    AABB box = make_AABB(start, end);
    AABB_clamp(box, parent);

    E_Int count = end - start;
    if (count <= Smesh::MAX_FACES_PER_BVH_LEAF) {
        return make_BVH_node(box, start, end, NULL, NULL);
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

    std::sort(bvh_fids.begin() + start, bvh_fids.begin() + end,
        [&] (E_Int fi, E_Int fj) {
            E_Float *fci = &fcenters[3*fi];
            E_Float *fcj = &fcenters[3*fj];
            return fci[dim] < fcj[dim];
        }
    );

    E_Int mid = start + count/2;

    BVH_node *left = make_BVH_subtree(start, mid, box);
    BVH_node *right = make_BVH_subtree(mid, end, box);

    return make_BVH_node(box, start, end, left, right);
}

void Smesh::make_BVH()
{
    bvh_fids.clear();
    bvh_fids.reserve(nf);
    for (E_Int i = 0; i < nf; i++) {
        bvh_fids.push_back(i);
    }

    bvh_root = make_BVH_subtree(0, nf, AABB_HUGE);
}

void Smesh::make_BVH(const std::set<E_Int> &fids)
{
    bvh_fids.clear();
    bvh_fids.reserve(fids.size());
    for (auto fid : fids) bvh_fids.push_back(fid);

    bvh_root = make_BVH_subtree(0, fids.size(), AABB_HUGE);
}

void Smesh::destroy_BVH(BVH_node *root)
{
    if (root == NULL) return;

    destroy_BVH(root->left);
    destroy_BVH(root->right);

    delete root;
}
