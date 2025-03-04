/*    
    Copyright 2013-2025 Onera.

    This file is part of Cassiopee.

    Cassiopee is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    Cassiopee is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Cassiopee.  If not, see <http://www.gnu.org/licenses/>.
*/
#pragma once

#include "common/common.h"

struct Mesh;

const E_Int normalIn_H[6] = {1, 0, 1, 0, 1, 0};

void H27_refine(E_Int hexa, Mesh *M);

void H27_reorder(E_Int hexa, Mesh *M);

void H18_refine(E_Int hexa, Mesh *M);

void H18_reorder(E_Int hexa, Mesh *M);

void reconstruct_quad(Mesh *M, E_Int hexa, E_Int *fids, E_Int crange, E_Int normalIn,
    E_Int NODE, E_Int pn[4]);

E_Int check_canon_hexa(E_Int hexa, Mesh *M);

void update_range_and_stride(Mesh *M, E_Int hexa, E_Int cpos, E_Int nchildren);
