#pragma once

#include <vector>

#include "xcore.h"

struct Hedge;
struct Vertex;

struct Cycle {
    Hedge *rep;
    E_Int inout;
    Vertex *left;
    Cycle *prev;
    Cycle *next;

    static E_Int INNER;
    static E_Int OUTER;
    static E_Int DEGEN;

    Cycle(Hedge *Rep);
    void print() const;

    static void set_inout(std::vector<Cycle *> &C);

    static void write_vertices(const char *fname, const std::vector<Cycle *> &C);
};