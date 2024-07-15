#pragma once

#include "../common/common.h"

extern const Int normalIn_Py[5];

struct Mesh;

void refine_pyra(Int pyra, Mesh *M);

void reorder_pyra(Int pyra, Mesh *M);

Int check_canon_pyra(Int pyra, Mesh *M);