#include "BVH.h"
#include "mesh.h"

// Creates the tree

static size_t leaf_count = 0;

BVH_node *BVH_create_node(AABB aabb, BVH_node *left_child, BVH_node *right_child,
    E_Int *ids, E_Int count, Mem_arena &arena)
{
    BVH_node *node = (BVH_node *)ARENA_ALLOC(sizeof(BVH_node));

    node->bbox = aabb;
    node->left = left_child;
    node->right = right_child;
    node->elements = (E_Int *)ARENA_ALLOC(count * sizeof(E_Int));
    assert(count > 0);
    memcpy(node->elements, ids, count * sizeof(E_Int));

    leaf_count++;

    return node;
}


BVH_node *BVH_create_node(AABB aabb, BVH_node *left_child,
    BVH_node *right_child, Mem_arena &arena)
{
    BVH_node *node = (BVH_node *)ARENA_ALLOC(sizeof(BVH_node));

    node->bbox = aabb;
    node->left = left_child;
    node->right = right_child;
    node->elements = NULL;

    return node;
}

BVH_node *BVH_create_node(E_Int *ids, E_Int count, const IMesh &M,
    Mem_arena &arena)
{
    AABB aabb(M, ids, count);

    if (count <= MAX_FACES_PER_LEAF) {
        return BVH_create_node(aabb, NULL, NULL, ids, count, arena);
    }

    E_Float dx = aabb.xmax - aabb.xmin;
    E_Float dy = aabb.ymax - aabb.ymin;
    E_Float dz = aabb.zmax - aabb.zmin;

    E_Float *pCoor;

    if (dx >= dy && dx >= dz) {
        pCoor = (E_Float *)M.X.data();
    } else if (dy >= dz) {
        assert(dy >= dx);
        pCoor = (E_Float *)M.Y.data();
    } else {
        assert(dz >= dx && dz >= dy);
        pCoor = (E_Float *)M.Z.data();
    }

    std::sort(ids, ids + count, [&] (E_Int i, E_Int j) 
    {
        E_Float cci = 0.0, ccj = 0.0;
        
        const auto &pni = M.F[i];
        for (const E_Int p : pni) cci += pCoor[p];
        cci /= pni.size();

        const auto &pnj = M.F[j];
        for (const E_Int p : pnj) ccj += pCoor[p];
        ccj /= pnj.size();
        
        return cci < ccj;
    });

    BVH_node *left = BVH_create_node(ids, count/2, M, arena);
    BVH_node *right = BVH_create_node(ids + count/2, count - count/2, M, arena);

    return BVH_create_node(aabb, left, right, arena);
}

BVH::BVH(const IMesh &M)
{
    if (M.skin.empty()) {
        puts("empty skin!");
        exit(1);
    }

    size_t NF = M.skin.size();

    printf("Number of faces: %zu\n", NF);

    arena.reserve(10 * NF * sizeof(E_Int));

    E_Int *ids = (E_Int *)ARENA_ALLOC(NF * sizeof(E_Int));

    for (size_t i = 0; i < NF; i++) ids[i] = M.skin[i];

    root = BVH_create_node(ids, NF, M, arena);

    printf("Leaves: %zu\n", leaf_count);

    arena.print_stats();
}