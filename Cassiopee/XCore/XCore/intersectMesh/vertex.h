#pragma once

#include <cstdio>

#include "xcore.h"

struct Hedge;

struct Vertex {
    E_Float x, y;
    Hedge *rep;
    E_Int id;
    Hedge *left;
    E_Int oid[2];

    Vertex(E_Float X, E_Float Y, E_Int Oid, E_Int color);

    Vertex(E_Float X, E_Float Y);

    inline void print() { printf("P" SF_D_ ": %f %f\n", id, x, y); }
};
