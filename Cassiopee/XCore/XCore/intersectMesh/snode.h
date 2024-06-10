#pragma once

#include <cstdio>

#include "segment.h"

struct Snode {
    Segment *key;
    void *inf;
    Snode *left;
    Snode *right;

    Snode(Segment *Key, void *Inf);

    inline void print()
    {
        printf("S" SF_D_ "(P" SF_D_ ", P" SF_D_ ")\n", key->id, key->p->id,
            key->q->id);
    }

    static void print_tree(Snode *root);
};
