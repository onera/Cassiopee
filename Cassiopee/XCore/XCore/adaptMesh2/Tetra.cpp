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

void refine_tetra(E_Int cell, AMesh *M)
{
  /*
  E_Int NODES[10], FACES[16];
  E_Int *BOT = FACES;
  E_Int *LFT = FACES + 4;
  E_Int *RGT = FACES + 8;
  E_Int *FRO = FACES + 12;
  E_Int *pf = get_facets(cell, M->nface, M->indPH);

  E_Int face, i0, reorient, *children, nchildren, local[6];
  Element **faceTree = M->faceTree;

  // BOT
  face = pf[0];
  reorient = get_reorient(face, cell, normalIn_T[0], M);
  children = faceTree[face]->children;
  for (E_Int i = 0; i < 4; i++) BOT[i] = children[i];
  NODES[0] = get_facets(face, M->ngon, M->indPG)[0];
  T6_get_ordered_data(M, NODES[0], reorient, BOT, local);
  assert(local[0] == NODES[0]);
  for (E_Int i = 1; i < 3; i++) NODES[i] = local[i];
  NODES[4] = local[3];
  NODES[5] = local[4];
  NODES[6] = local[5];

  // LFT
  face = pf[1];
  reorient = get_reorient(face, cell, normalIn_T[1], M);
  children = faceTree[face]->children;
  for (E_Int i = 0; i < 4; i++) LFT[i] = children[i]; 
  T6_get_ordered_data(M, NODES[0], reorient, LFT, local);
  assert(local[5] == NODES[6]);
  NODES[3] = local[1];
  NODES[7] = local[4];
  NODES[8] = local[3];

  // RGT
  face = pf[2];
  reorient = get_reorient(face, cell, normalIn_T[2], M);
  children = faceTree[face]->children;
  for (E_Int i = 0; i < 4; i++) RGT[i] = children[i]; 
  T6_get_ordered_data(M, NODES[1], reorient, RGT, local);
  assert(local[0] == NODES[1]);
  assert(local[1] == NODES[3]); 
  assert(local[2] == NODES[2]);
  assert(local[4] == NODES[7]);
  assert(local[5] == NODES[5]);
  NODES[9] = local[3];

  // FRO
  face = pf[3];
  reorient = get_reorient(face, cell, normalIn_T[3], M);
  children = faceTree[face]->children;
  for (E_Int i = 0; i < 4; i++) FRO[i] = children[i]; 
  T6_get_ordered_data(M, NODES[0], reorient, FRO, local);
  assert(local[0] == NODES[0]);
  assert(local[1] == NODES[1]);
  assert(local[2] == NODES[3]);
  assert(local[3] == NODES[4]);
  assert(local[4] == NODES[9]);
  assert(local[5] == NODES[8]);

  // Set internal faces in ngon
  E_Int *ptr = &M->indPG[M->nfaces];

  // nfaces
  ptr[1] = ptr[0] + 3;
  M->ngon[*ptr  ] = NODES[4];  M->ngon[*ptr+1] = NODES[8];
  M->ngon[*ptr+2] = NODES[6];
  ptr++;
  
  // nfaces+1
  ptr[1] = ptr[0] + 3;
  M->ngon[*ptr  ] = NODES[4];  M->ngon[*ptr+1] = NODES[9];
  M->ngon[*ptr+2] = NODES[5];
  ptr++; 

  // nfaces+2
  ptr[1] = ptr[0] + 3;
  M->ngon[*ptr  ] = NODES[6];  M->ngon[*ptr+1] = NODES[5];
  M->ngon[*ptr+2] = NODES[7];
  ptr++;
  
  // nfaces+3
  ptr[1] = ptr[0] + 3;
  M->ngon[*ptr  ] = NODES[8];  M->ngon[*ptr+1] = NODES[9];
  M->ngon[*ptr+2] = NODES[7];
  ptr++;

  // Connect nodes 4 and 7 -> four new faces
  // TODO(Imad): test the other two choices of diagonal
  // nfaces+4 (const)
  ptr[1] = ptr[0] + 3; 
  M->ngon[*ptr  ] = NODES[6];  M->ngon[*ptr+1] = NODES[4];
  M->ngon[*ptr+2] = NODES[7];
  ptr++;
  
  // nfaces+5
  ptr[1] = ptr[0] + 3;
  M->ngon[*ptr  ] = NODES[4];  M->ngon[*ptr+1] = NODES[9];
  M->ngon[*ptr+2] = NODES[7];
  ptr++;
  
  // nfaces+6 (const)
  ptr[1] = ptr[0] + 3;
  M->ngon[*ptr  ] = NODES[4];  M->ngon[*ptr+1] = NODES[8];
  M->ngon[*ptr+2] = NODES[7];
  ptr++;
  
  // nfaces+7
  ptr[1] = ptr[0] + 3;
  M->ngon[*ptr  ] = NODES[4];  M->ngon[*ptr+1] = NODES[7];
  M->ngon[*ptr+2] = NODES[5];

  // Assemble children
  
  // First child replaces cell
  ptr = &M->indPH[cell];
  M->nface[*ptr  ] = BOT[0];      M->nface[*ptr+1] = LFT[0];
  M->nface[*ptr+2] = M->nfaces;   M->nface[*ptr+3] = FRO[0];

  ptr = &M->indPH[M->ncells]; 

  // ncells
  ptr[1] = ptr[0] + 4;
  M->nface[*ptr  ] = BOT[1];      M->nface[*ptr+1] = M->nfaces+1;
  M->nface[*ptr+2] = RGT[0];      M->nface[*ptr+3] = FRO[1];
  ptr++;

  // ncells+1
  ptr[1] = ptr[0] + 4;
  M->nface[*ptr  ] = BOT[2];      M->nface[*ptr+1] = LFT[2];
  M->nface[*ptr+2] = RGT[2];      M->nface[*ptr+3] = M->nfaces+2;
  ptr++;
  
  // ncells+2
  ptr[1] = ptr[0] + 4;
  M->nface[*ptr  ] = M->nfaces+3; M->nface[*ptr+1] = LFT[1];
  M->nface[*ptr+2] = RGT[1];      M->nface[*ptr+3] = FRO[2];
  ptr++;

  // Octahedron -> 4 new tetra (clockwise)
  // ncells+3
  ptr[1] = ptr[0] + 4;
  M->nface[*ptr  ] = M->nfaces+4; M->nface[*ptr+1] = LFT[3];
  M->nface[*ptr+2] = M->nfaces+6; M->nface[*ptr+3] = M->nfaces;
  ptr++;
 
  // ncells+4
  ptr[1] = ptr[0] + 4;
  M->nface[*ptr  ] = M->nfaces+5; M->nface[*ptr+1] = M->nfaces+6;
  M->nface[*ptr+2] = M->nfaces+3; M->nface[*ptr+3] = FRO[3];
  ptr++;

  // ncells+5
  ptr[1] = ptr[0] + 4;
  M->nface[*ptr  ] = M->nfaces+2; M->nface[*ptr+1] = M->nfaces+4;
  M->nface[*ptr+2] = M->nfaces+7; M->nface[*ptr+3] = BOT[3];
  ptr++;
  
  // ncells+6
  ptr[1] = ptr[0] + 4;
  M->nface[*ptr  ] = M->nfaces+1; M->nface[*ptr+1] = M->nfaces+7;
  M->nface[*ptr+2] = RGT[3]; M->nface[*ptr+3] = M->nfaces+5;

  // Set new cells in tree

  E_Int next_level = M->cellTree[cell]->level + 1;

  // Parent
  Element *Elem = M->cellTree[cell];
  Elem->children = (E_Int *)XMALLOC(8 * sizeof(E_Int));
  Elem->children[0] = cell;
  for (E_Int i = 1; i < 8; i++) Elem->children[i] = M->ncells+i-1;
  Elem->nchildren = 8;
  Elem->parent = cell;
  Elem->position = 0;
  Elem->type = TETRA;
  Elem->level = next_level;
  Elem->state = 1;

  // Tetra children
  for (E_Int i = 0; i < 7; i++) {
    Element *Elem = M->cellTree[M->ncells+i];
    Elem->children = NULL;
    Elem->nchildren = 0;
    Elem->parent = cell;
    Elem->position = i+1;
    Elem->type = TETRA;
    Elem->level = next_level;
    Elem->state = 0;
  }

  //E_Int next_level = M->cellTree[cell]->level + 1;
  //set_parent_elem(M->cellTree, cell, 8, TETRA, M->ncells);
  //for (E_Int i = 0; i < 7; i++)
  //  set_child_elem(M->cellTree, cell, i, TETRA, next_level, M->ncells);

  // Set external faces owns and neis
  update_external_own_nei(cell, M);

  // Set internal faces in face tree
  for (E_Int i = 0; i < 8; i++)
    set_new_face_in_tree(faceTree, M->nfaces+i, TRI, next_level);

  // Set owns and neis of internal faces
  M->owner[M->nfaces]   = M->ncells+3; M->neigh[M->nfaces]   = cell;
  M->owner[M->nfaces+1] = M->ncells+0; M->neigh[M->nfaces+1] = M->ncells+6;
  M->owner[M->nfaces+2] = M->ncells+1; M->neigh[M->nfaces+2] = M->ncells+5;
  M->owner[M->nfaces+3] = M->ncells+4; M->neigh[M->nfaces+3] = M->ncells+2;

  M->owner[M->nfaces+4] = M->ncells+5; M->neigh[M->nfaces+4] = M->ncells+3;
  M->owner[M->nfaces+5] = M->ncells+6; M->neigh[M->nfaces+5] = M->ncells+4;
  M->owner[M->nfaces+6] = M->ncells+4; M->neigh[M->nfaces+6] = M->ncells+3;
  M->owner[M->nfaces+7] = M->ncells+6; M->neigh[M->nfaces+7] = M->ncells+5;

  check_canon_tetra(cell, M);
  for (E_Int i = 0; i < 7; i++)
    check_canon_tetra(M->ncells+i, M);

  M->ncells += 7;
  M->nfaces += 8;
  */
}

void reorder_tetra(E_Int i, E_Int nf, E_Int *pf, AMesh *M)
{
  E_Int common[3], map[3];
  E_Int bot = pf[0];
  E_Int *pn = get_facets(bot, M->ngon, M->indPG);

  for (E_Int j = 0; j < 3; j++) map[j] = pn[j];
  E_Int reorient = get_reorient(bot, i, normalIn_T[0], M);
  if (reorient) std::swap(map[1], map[2]);

  E_Int lft, rgt, fro;
  lft = rgt = fro = -1;

  for (E_Int j = 1; j < 4; j++) {
    E_Int face = pf[j];

    for (E_Int k = 0; k < 3; k++) common[k] = 0;

    E_Int *pnn = get_facets(face, M->ngon, M->indPG);

    for (E_Int k = 0; k < 3; k++) {
      E_Int point = pnn[k];

      for (E_Int l = 0; l < 3; l++) {
        if (map[l] == point) {
          common[l] = 1;
          break;
        }
      }
    }

    if      (common[0] && common[2]) lft = face;
    else if (common[1] && common[2]) rgt = face;
    else if (common[1] && common[0]) fro = face;
    else assert(0);
  }
  assert(lft != -1 && rgt != -1 && fro != -1);

  pf[1] = lft;
  pf[2] = rgt;
  pf[3] = fro;
}

E_Int check_canon_tetra(E_Int cell, AMesh *M)
{
  E_Int NODES[4] = {-1, -1, -1, -1};

  E_Int *pf = get_facets(cell, M->nface, M->indPH);

  E_Int face, i0, reorient, *pn, local[3];

  // BOT (In)
  face = pf[0];
  pn = get_facets(face, M->ngon, M->indPG);
  i0 = 0;
  reorient = get_reorient(face, cell, normalIn_T[0], M);
  reorder_tri(local, pn, reorient, i0);
  for (E_Int i = 0; i < 3; i++) NODES[i] = local[i];

  // LFT (Out)
  face = pf[1];
  pn = get_facets(face, M->ngon, M->indPG);
  i0 = Get_pos(NODES[0], pn, 3);
  reorient = get_reorient(face, cell, normalIn_T[1], M);
  reorder_tri(local, pn, reorient, i0);
  assert(local[0] == NODES[0]);
  assert(local[2] == NODES[2]);
  NODES[3] = local[1];

  // RGT (In)
  face = pf[2];
  pn = get_facets(face, M->ngon, M->indPG);
  i0 = Get_pos(NODES[1], pn, 3);
  reorient = get_reorient(face, cell, normalIn_T[2], M);
  reorder_tri(local, pn, reorient, i0);
  assert(local[0] == NODES[1]);
  assert(local[1] == NODES[3]);
  assert(local[2] == NODES[2]);

  // FRO (Out)
  face = pf[3];
  pn = get_facets(face, M->ngon, M->indPG);
  i0 = Get_pos(NODES[0], pn, 3);
  reorient = get_reorient(face, cell, normalIn_T[3], M);
  reorder_tri(local, pn, reorient, i0);
  assert(local[0] == NODES[0]);
  assert(local[1] == NODES[1]);
  assert(local[2] == NODES[3]);

  return 1;
}

void make_ref_data_tetra(E_Int cell, AMesh *M, E_Float *pM, const pDirs &Dirs)
{}

void make_pdirs_tetra(E_Int cell, AMesh *M, pDirs &Dirs)
{}