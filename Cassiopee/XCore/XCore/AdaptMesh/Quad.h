/*    
    Copyright 2013-2024 Onera.

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
#include "Mesh.h"

Int Q9_refine(Int quad, Mesh *M);

Int Q6_refine(Int quad, Mesh *M);

void Q4_reorder(Int *pn, Int reorient, Int i0, Int local[4]);

inline
void refine_face_iso(Int face, Mesh *M)
{
    switch (M->ftype[face]) {
        case QUAD:
            Q9_refine(face, M);
            break;
        default:
            assert(0);
            break;
    }
}