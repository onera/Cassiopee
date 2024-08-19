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

static inline
void refine_cell_iso(E_Int cell, Mesh *M)
{
    switch (M->ctype[cell]) {
        case HEXA:
            H27_refine(cell, M);
            break;
        default:
            assert(0);
            break;
    }
}

void Mesh_refine_iso(Mesh *M, std::vector<E_Int> &ref_cells,
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
            ref_faces.pop_back();
            refine_face_iso(face, M);
        }

        while (!ref_cells.empty() && M->clevel[ref_cells.back()] == level) {
            E_Int cell = ref_cells.back();
            ref_cells.pop_back();
            refine_cell_iso(cell, M);
        }
    }
}

/*
void Mesh_refine_iso(Mesh *M, std::vector<E_Int> &ref_cells,
    std::vector<E_Int> &ref_faces, std::set<UEdge> &ref_edges)
{
    // Adapt the cells that do not need face adaptation

    if (ref_faces.empty()) {
        for (E_Int cell : ref_cells) refine_cell_iso(cell, M);
        
        return;
    }

    // Adapt the faces that are not linked to local cells

    if (ref_cells.empty()) {
        for (E_Int face : ref_faces) refine_face_iso(face, M);

        return;
    }

    // Refine first cells whose level is less than the smallest ref_face level
    // Note(Imad): the faces are pre-sorted by increasing level

    E_Int min_flvl = M->flevel[ref_faces[0]];

    // Refine the "lagging" cells

    E_Int cell_start = 0;

    while (M->clevel[ref_cells[cell_start]] < min_flvl) {
        refine_cell_iso(ref_cells[cell_start], M);
        cell_start++;
    }

    // Now refine lagging faces
    
    E_Int min_clvl = M->clevel[ref_cells[cell_start]];

    E_Int face_start = 0;

    while (M->flevel[ref_faces[face_start]] < min_clvl) {
        refine_face_iso(ref_faces[face_start], M);
        face_start++;
    }


    // At this point the remaining ref faces/cells should be at the same level

    E_Int cells_left = (E_Int)ref_cells.size() - cell_start;
    E_Int faces_left = (E_Int)ref_faces.size() - face_start;

    // Refine faces and cells, by increasing level

    E_Int current_lvl = M->clevel[ref_cells[cell_start]];

    while (cells_left || faces_left) {

        while (faces_left) {
            
            // Break early if refinement would overstep current level
            
            if (M->flevel[ref_faces[face_start]] > current_lvl) break;
            
            E_Int face = ref_faces[face_start];

            refine_face_iso(face, M);

            face_start++;
            faces_left--;
        }

        while (cells_left) {

            // Break early if refinement would overstep current level
            
            if (M->clevel[ref_cells[cell_start]] > current_lvl) break;

            E_Int cell = ref_cells[cell_start];

            refine_cell_iso(cell, M);

            cell_start++;
            cells_left--;
        }

        current_lvl++;
    }
}
*/