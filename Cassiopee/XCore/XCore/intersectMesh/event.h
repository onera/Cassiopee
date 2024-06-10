#pragma once

#include <vector>

#include "xcore.h"

struct Vertex;
struct Segment;

struct Event {
    Vertex *key;
    Segment *inf;
    Event *left;
    Event *right;

    Event(E_Float x, E_Float y, E_Int oid, E_Int color);
    Event(E_Float x, E_Float y);

    void inorder(std::vector<Vertex *> &V) const;
};
