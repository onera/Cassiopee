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
void refine_cell_dir(E_Int cell, Mesh *M)
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
void refine_face_dir(E_Int face, E_Int pattern, Mesh *M)
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
E_Int get_face_pattern(E_Int fid, Mesh *M)
{
    E_Int own = M->owner[fid];
    E_Int fpos = -1;
    E_Int *cell = Mesh_get_cell(M, own);
    E_Int *crange = Mesh_get_crange(M, own);

    for (E_Int j = 0; j < M->cstride[own] && fpos == -1; j++) {
        E_Int *pf = cell + 4*j;
        for (E_Int k = 0; k < crange[j]; k++) {
            E_Int face = pf[k];
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

void Mesh_refine_dir(Mesh *M, std::vector<E_Int> &ref_cells,
    std::vector<E_Int> &ref_faces, std::set<UEdge> &ref_edges)
{
    std::set<E_Int> levelset;

    for (E_Int cell : ref_cells) {
        levelset.insert(M->clevel[cell]);
    }

    for (E_Int face : ref_faces) {
        levelset.insert(M->flevel[face]);
    }

    std::vector<E_Int> levels;
    for (E_Int level : levelset) levels.push_back(level);
    std::sort(levels.begin(), levels.end());

    std::reverse(ref_cells.begin(), ref_cells.end());
    std::reverse(ref_faces.begin(), ref_faces.end());

    for (E_Int level : levels) {
        while (!ref_faces.empty() && M->flevel[ref_faces.back()] == level) {
            E_Int face = ref_faces.back();
            E_Int pattern = get_face_pattern(face, M);
            ref_faces.pop_back();
            refine_face_dir(face, pattern, M);
        }

        while (!ref_cells.empty() && M->clevel[ref_cells.back()] == level) {
            E_Int cell = ref_cells.back();
            ref_cells.pop_back();
            refine_cell_dir(cell, M);
        }
    }
}