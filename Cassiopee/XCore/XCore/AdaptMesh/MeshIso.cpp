#include "Mesh.h"
#include "common/mem.h"
#include "Hexa.h"
#include "Quad.h"
#include "Edge.h"

static inline
void refine_cell_iso(Int cell, Mesh *M)
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

void Mesh_refine_iso(Mesh *M, std::vector<Int> &ref_cells,
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
            ref_faces.pop_back();
            refine_face_iso(face, M);
        }

        while (!ref_cells.empty() && M->clevel[ref_cells.back()] == level) {
            Int cell = ref_cells.back();
            ref_cells.pop_back();
            refine_cell_iso(cell, M);
        }
    }
}

/*
void Mesh_refine_iso(Mesh *M, std::vector<Int> &ref_cells,
    std::vector<Int> &ref_faces, std::set<UEdge> &ref_edges)
{
    // Adapt the cells that do not need face adaptation

    if (ref_faces.empty()) {
        for (Int cell : ref_cells) refine_cell_iso(cell, M);
        
        return;
    }

    // Adapt the faces that are not linked to local cells

    if (ref_cells.empty()) {
        for (Int face : ref_faces) refine_face_iso(face, M);

        return;
    }

    // Refine first cells whose level is less than the smallest ref_face level
    // Note(Imad): the faces are pre-sorted by increasing level

    Int min_flvl = M->flevel[ref_faces[0]];

    // Refine the "lagging" cells

    Int cell_start = 0;

    while (M->clevel[ref_cells[cell_start]] < min_flvl) {
        refine_cell_iso(ref_cells[cell_start], M);
        cell_start++;
    }

    // Now refine lagging faces
    
    Int min_clvl = M->clevel[ref_cells[cell_start]];

    Int face_start = 0;

    while (M->flevel[ref_faces[face_start]] < min_clvl) {
        refine_face_iso(ref_faces[face_start], M);
        face_start++;
    }


    // At this point the remaining ref faces/cells should be at the same level

    Int cells_left = (Int)ref_cells.size() - cell_start;
    Int faces_left = (Int)ref_faces.size() - face_start;

    // Refine faces and cells, by increasing level

    Int current_lvl = M->clevel[ref_cells[cell_start]];

    while (cells_left || faces_left) {

        while (faces_left) {
            
            // Break early if refinement would overstep current level
            
            if (M->flevel[ref_faces[face_start]] > current_lvl) break;
            
            Int face = ref_faces[face_start];

            refine_face_iso(face, M);

            face_start++;
            faces_left--;
        }

        while (cells_left) {

            // Break early if refinement would overstep current level
            
            if (M->clevel[ref_cells[cell_start]] > current_lvl) break;

            Int cell = ref_cells[cell_start];

            refine_cell_iso(cell, M);

            cell_start++;
            cells_left--;
        }

        current_lvl++;
    }
}
*/