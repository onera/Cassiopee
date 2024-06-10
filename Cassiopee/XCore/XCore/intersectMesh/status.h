#pragma once

#include <cstdio>

#include "snode.h"

struct Status {
    Snode *root;
    E_Float rx, ry;

    Status();

    Snode *insert(Segment *key, void *inf = NULL);

    Snode *lookup(Segment *seg);

    Snode *locate(Segment *seg);

    Snode *pred(Snode *sit);
    
    Snode *pred(Segment *seg);
    
    Snode *succ(Snode *sit);
    
    Snode *succ(Segment *seg);

    void erase(Segment *seg);

    void print();
    
    Snode *insert_(Snode *& root, Segment *key, void *inf);
    
    Snode *lookup_(Snode *root, Segment *seg);

    Snode *locate_(Snode *root, Segment *seg);

    void pred_(Snode *root, Segment *seg, Snode *&pre);
    
    void succ_(Snode *root, Segment *seg, Snode *&suc);

    Snode *erase_(Snode *root, Segment *seg);
};
