#pragma once

#include <cstddef>
#include <vector>

#include "xcore.h"

struct Vertex;
struct Segment;
struct Event;

struct Queue {
    Event *root;

    E_Int nelem;

    Queue();

    Event *insert(E_Float X, E_Float Y, E_Int oid, E_Int color);

    Event *insert(E_Float X, E_Float Y);

    Event *lookup(Vertex *key);
    
    Event *lookup(E_Float x, E_Float y);
    
    inline bool empty() { return root == NULL; }

    Event *min();

    void erase(Event *event);

    void erase(Vertex *p);

    void inorder(std::vector<Vertex *> &V) const;
    
    Event *insert_(Event *&root, E_Float x, E_Float y, E_Int oid, E_Int color);

    Event *insert_(Event *&root, E_Float x, E_Float y);

    Event *lookup_(Event *root, E_Float x, E_Float y);

    Event *erase_(Event *root, Vertex *p);
};
