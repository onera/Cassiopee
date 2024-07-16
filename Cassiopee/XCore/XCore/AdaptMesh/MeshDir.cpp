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
#include "Mesh.h"
#include "common/mem.h"
#include "Hexa.h"
#include "Quad.h"
#include "Edge.h"

#define ISO 0
#define DIR 1

static inline
void refine_cell_dir(Int cell, Mesh *M)
{
    switch (M->ctype[cell]) {
        case HEXA:
            H18_refine(cell, M);
            break;
        default:
            assert(0);
            break;
    }
}

static inline
void refine_face_dir(Int face, Int pattern, Mesh *M)
{
    switch (M->ftype[face]) {
        case QUAD: {
            if (pattern == ISO) Q9_refine(face, M);
            else Q6_refine(face, M);
            break;
        }
        default:
            assert(0);
            break;
    }
}

static inline
Int get_face_pattern(Int fid, Mesh *M)
{
    Int own = M->owner[fid];
    Int cref = M->cref[own];
    Int fpos = -1;
    Int *cell = Mesh_get_cell(M, own);
    Int *crange = Mesh_get_crange(M, own);
    for (Int j = 0; j < M->cstride[own] && fpos == -1; j++) {
        Int *pf = cell + 4*j;
        for (Int k = 0; k < crange[j]; k++) {
            Int face = pf[k];
            if (face == fid) {
                fpos = j;
                assert(k == 0);
                break;
            }
        }
    }

    assert(fpos != -1);

    if (fpos == 0 || fpos == 1) return ISO;
    return DIR;
}

void Mesh_refine_dir(Mesh *M, std::vector<Int> &ref_cells,
    std::vector<Int> &ref_faces, std::set<UEdge> &ref_edges)
{
    std::set<Int> levelset;

    for (Int cell : ref_cells) {
        levelset.insert(M->clevel[cell]);
    }

    for (Int face : ref_faces) {
        levelset.insert(M->flevel[face]);
    }

    std::vector<Int> levels;
    for (Int level : levelset) levels.push_back(level);
    std::sort(levels.begin(), levels.end());

    std::reverse(ref_cells.begin(), ref_cells.end());
    std::reverse(ref_faces.begin(), ref_faces.end());

    for (Int level : levels) {
        while (!ref_faces.empty() && M->flevel[ref_faces.back()] == level) {
            Int face = ref_faces.back();
            Int pattern = get_face_pattern(face, M);
            ref_faces.pop_back();
            refine_face_dir(face, pattern, M);
        }

        while (!ref_cells.empty() && M->clevel[ref_cells.back()] == level) {
            Int cell = ref_cells.back();
            ref_cells.pop_back();
            refine_cell_dir(cell, M);
        }
    }
}