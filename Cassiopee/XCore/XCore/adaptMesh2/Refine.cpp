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
#include "Proto.h"
#include <set>

void resize_data_for_refinement(AMesh *M, size_t nref_cells, size_t nref_faces)
{
  E_Int cell_incr = nref_cells * 7;
  E_Int face_incr = nref_faces * 3 // OK for quads and tris
                  + nref_cells * 12; // max 13 internal faces per refined cell
  E_Int nnew_cells = M->ncells + cell_incr;
  E_Int nnew_faces = M->nfaces + face_incr;

  // max 5 new points per refined quad + nref_cells centroids
  E_Int nnew_points = M->npoints + nref_faces*5 + nref_cells;
 
  M->cellTree->resize(nnew_cells);
  M->faceTree->resize(nnew_faces);
  
  M->ngon  = (E_Int *)  XRESIZE(M->ngon,  (4*nnew_faces) * sizeof(E_Int));
  M->indPG = (E_Int *)  XRESIZE(M->indPG, (nnew_faces+1) * sizeof(E_Int));
  M->nface = (E_Int *)  XRESIZE(M->nface, (6*nnew_cells) * sizeof(E_Int));
  M->indPH = (E_Int *)  XRESIZE(M->indPH, (nnew_cells+1) * sizeof(E_Int));
  M->owner = (E_Int *)  XRESIZE(M->owner, (nnew_faces)   * sizeof(E_Int));
  M->neigh = (E_Int *)  XRESIZE(M->neigh, (nnew_faces)   * sizeof(E_Int));
  M->x     = (E_Float *)XRESIZE(M->x,     (nnew_points)  * sizeof(E_Float));
  M->y     = (E_Float *)XRESIZE(M->y,     (nnew_points)  * sizeof(E_Float));
  M->z     = (E_Float *)XRESIZE(M->z,     (nnew_points)  * sizeof(E_Float));  

  for (E_Int i = M->nfaces; i < nnew_faces; i++) {
    M->owner[i] = -1;
    M->neigh[i] = -1;
  }
}

void update_external_pe_after_ref(E_Int cell, AMesh *M)
{
  Children *children = M->cellTree->children(cell);

  for (E_Int i = 0; i < children->n; i++) {
    E_Int child = children->pc[i];

    for (E_Int j = M->indPH[child]; j < M->indPH[child+1]; j++) {
      E_Int face = M->nface[j];
      if (M->owner[face] == cell) {
        M->owner[face] = child;
      } else if (M->neigh[face] == cell) {
        M->neigh[face] = child;
      }
    }
  }
}

void get_ref_faces_and_cells(AMesh *M, std::vector<E_Int> &ref_faces,
  std::vector<E_Int> &ref_cells)
{
  ref_cells.clear();
  ref_faces.clear();

  for (E_Int i = 0; i < M->ncells; i++) {
    if (M->ref_data[i] > 0) {
      ref_cells.push_back(i);
    }
  }

  // Analyse faces

  // Boundary faces
  for (E_Int i = 0; i < M->nbc; i++) {
    E_Int *ptlist = M->ptlists[i];
    
    for (E_Int j = 0; j < M->bcsizes[i]; j++) {
      E_Int face = ptlist[j];
      
      E_Int own = M->owner[face];

      E_Int flvl = M->faceTree->level(face);
      E_Int olvl = M->cellTree->level(own);
      assert(olvl == flvl);
      E_Int oval = olvl + M->ref_data[own];

      if (oval > flvl) {
        ref_faces.push_back(face);
      }
    }
  }

  // Proc faces
  for (E_Int i = 0; i < M->npatches; i++) {
    Patch *P = &M->patches[i];

    P->sbuf_i = (E_Int *)XRESIZE(P->sbuf_i, P->nf * sizeof(E_Int));
    P->rbuf_i = (E_Int *)XRESIZE(P->rbuf_i, P->nf * sizeof(E_Int));

    for (E_Int j = 0; j < P->nf; j++) {
      E_Int face = P->pf[j];
      E_Int own = M->owner[face];
      P->sbuf_i[j] = M->ref_data[own] + M->cellTree->level(own);
    }

    MPI_Isend(P->sbuf_i, P->nf, MPI_INT, P->nei, M->pid, MPI_COMM_WORLD,
      &M->req[M->nrq++]);
    MPI_Irecv(P->rbuf_i, P->nf, MPI_INT, P->nei, P->nei, MPI_COMM_WORLD,
      &M->req[M->nrq++]);
  }

  Comm_waitall(M);

  for (E_Int i = 0; i < M->npatches; i++) {
    Patch *P = &M->patches[i];

    for (E_Int j = 0; j < P->nf; j++) {
      E_Int face = P->pf[j];
      E_Int own = M->owner[face];

      E_Int flvl = M->faceTree->level(face);
      E_Int oval = M->ref_data[own] + M->cellTree->level(own);
      E_Int nval = P->rbuf_i[j];

      if (oval > flvl || nval > flvl) {
        ref_faces.push_back(face);
      }
    }
  }

  // Internal faces
  for (E_Int i = 0; i < M->nfaces; i++) {
    E_Int nei = M->neigh[i];

    if (nei == -1) continue;
    
    E_Int flvl = M->faceTree->level(i);
    
    E_Int own = M->owner[i];
    
    E_Int lo = M->cellTree->level(own);
    E_Int ro = M->ref_data[own];
    E_Int oval = lo + ro;
    
    E_Int ln = M->cellTree->level(nei);
    E_Int rn = M->ref_data[nei];
    E_Int nval = ln + rn;
    
    if (oval > flvl || nval > flvl) {
      ref_faces.push_back(i);
      continue;
    }
  }
}

static
void assign_pattern_to_ref_face(E_Int face, E_Int &pattern, AMesh *M)
{
  if (!M->mode_2D) return;

  E_Int own = M->owner[face];
  E_Int nf = -1;
  E_Int *pf = get_cell(own, nf, M->nface, M->indPH);

  E_Int fpos = -1;
    
  E_Int flvl = M->faceTree->level(face);
  E_Int olvl = M->cellTree->level(own);

  if (flvl == olvl)
    fpos = Get_pos(face, pf, nf);
  else
    fpos = Get_pos(master_face(face, M), pf, nf);
    
  assert(fpos != -1);

  // BOT and TOP: ISO
  // Other faces: DIR

  if (fpos == 0 || fpos == 1) pattern = ISO;
  else pattern = DIR;
}

static
void refine_cell(E_Int cell, AMesh *M)
{
  if (M->mode_2D) H18_refine(cell, M);
  else H27_refine(cell, M);
}

static
void refine_face(E_Int face, E_Int pattern, AMesh *M)
{
  if (!M->mode_2D) {
    Q9_refine(face, M);
  } else {
    if      (pattern == DIR) Q6_refine(face, M);
    else if (pattern == ISO) Q9_refine(face, M);
    else assert(0);
  }
}

void refine_mesh(AMesh *M, const std::vector<E_Int> &ref_faces,
  const std::vector<E_Int> &ref_cells)
{
  // First update the cells that do not need face refinement
  if (ref_faces.size() == 0) {
    for (auto cell: ref_cells) refine_cell(cell, M);

    update_boundary_faces(M);

    // TODO(Imad): get rid of this resize
    if (M->ncells < M->cellTree->nelem_)
      M->cellTree->resize(M->ncells);
    if (M->nfaces < M->faceTree->nelem_)
      M->faceTree->resize(M->nfaces);

    return;
  }

  if (ref_cells.size() == 0) {
    for (auto face: ref_faces) {
      E_Int pattern;
      assign_pattern_to_ref_face(face, pattern, M);
      refine_face(face, pattern, M);
    }

    update_boundary_faces(M);

    // TODO(Imad): get rid of this resize
    if (M->ncells < M->cellTree->nelem_)
      M->cellTree->resize(M->ncells);
    if (M->nfaces < M->faceTree->nelem_)
      M->faceTree->resize(M->nfaces);

    return;
  }

  // Ensure cells with a smaller level than the smallest ref_face level
  // get refined first

  E_Int min_flvl = M->faceTree->level(ref_faces[0]);

  E_Int cell_start = 0; // Eliminate 'lagging' cells

  while (M->cellTree->level(ref_cells[cell_start]) < min_flvl) {
    refine_cell(ref_cells[cell_start], M);
    cell_start += 1;
  }

  // Now refine lagging faces

  E_Int min_clvl = M->cellTree->level(ref_cells[cell_start]);

  E_Int face_start = 0;

  while (M->faceTree->level(ref_faces[face_start]) < min_clvl) {
    E_Int pattern;
    E_Int face = ref_faces[face_start];
    assign_pattern_to_ref_face(face, pattern, M);
    refine_face(face, pattern, M);
    face_start += 1;
  }
 
  E_Int cells_left = ref_cells.size() - cell_start;
  E_Int faces_left = ref_faces.size() - face_start;

  // Now ref_faces and ref_cells should be at the same level

  E_Int cur_lvl = M->cellTree->level(ref_cells[cell_start]);

  while (cells_left || faces_left) {

    while (faces_left) {
      if (M->faceTree->level(ref_faces[face_start]) > cur_lvl)
        break;
      
      E_Int pattern;
      E_Int face = ref_faces[face_start];
      assign_pattern_to_ref_face(face, pattern, M);
      refine_face(face, pattern, M);
      face_start += 1;
      faces_left -= 1;
    }

    while (cells_left) {
      if (M->cellTree->level(ref_cells[cell_start]) > cur_lvl)
        break;

      E_Int cell = ref_cells[cell_start];
      refine_cell(cell, M);
      cell_start += 1;
      cells_left -= 1;
    }

    cur_lvl += 1;
  }

  update_boundary_faces(M);

  // TODO(Imad): get rid of this resize
  if (M->ncells < M->cellTree->nelem_)
    M->cellTree->resize(M->ncells);
  if (M->nfaces < M->faceTree->nelem_)
    M->faceTree->resize(M->nfaces);

}