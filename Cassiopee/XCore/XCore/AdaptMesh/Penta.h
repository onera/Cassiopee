#pragma once

#include "../common/common.h"

extern const Int normalIn_Pe[5];

struct Mesh;

void refine_penta(Int penta, Mesh *M);

void reorder_penta(Int penta, Mesh *M);

Int check_canon_penta(Int penta, Mesh *M);