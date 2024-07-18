#pragma once

#include "common/common.h"

extern const Int normalIn_T[4];

struct Mesh;

void refine_tetra(Int tetra, Mesh *M);

void reorder_tetra(Int tetra, Mesh *M);

Int check_canon_tetra(Int tetra, Mesh *M);