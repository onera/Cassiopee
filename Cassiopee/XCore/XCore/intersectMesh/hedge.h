#pragma once

#include <vector>

#include "xcore.h"

struct Vertex;
struct Face;
struct Cycle;

struct Hedge {
    Vertex *orig;
    Hedge *twin;
    Hedge *prev;
    Hedge *next;
    Face *left;
    E_Int color;
    Cycle *cycle;

    Hedge(Vertex *Orig);

    static E_Int cmp_cwise(const Hedge *h, const Hedge *w);
    static void sort_cwise(std::vector<Hedge *> &H, E_Int start, E_Int end);
    static void sort_ccwise(std::vector<Hedge *> &H, E_Int start, E_Int end);
};