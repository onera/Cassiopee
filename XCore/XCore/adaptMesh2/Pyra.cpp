#include "Proto.h"

void refine_pyra(E_Int cell, AMesh *M)
{
  /*
  E_Int NODES[14], FACES[20];
  for (E_Int i = 0; i < 14; i++) NODES[i] = -1;
  E_Int *BOT = FACES;
  E_Int *LFT = FACES + 4;
  E_Int *RGT = FACES + 8;
  E_Int *FRO = FACES + 12;
  E_Int *BCK = FACES + 16;
  E_Int *pf = get_facets(cell, M->nface, M->indPH);

  E_Int face, i0, reorient, *children, nchildren, local[9];
  Element **faceTree = M->faceTree;

  // BOT (In)
  face = pf[0];
  reorient = get_reorient(face, cell, normalIn_Py[0], M);
  children = faceTree[face]->children;
  for (E_Int i = 0; i < 4; i++) BOT[i] = children[i];
  NODES[0] = get_facets(face, M->ngon, M->indPG)[0];
  Q9_get_ordered_data(M, NODES[0], reorient, BOT, local);
  assert(NODES[0] == local[0]);
  for (E_Int i = 1; i < 4; i++) NODES[i] = local[i];
  NODES[5] = local[4];
  NODES[6] = local[5];
  NODES[7] = local[6];
  NODES[8] = local[7];
  NODES[9] = local[8];

  // LFT (In)
  face = pf[1];
  reorient = get_reorient(face, cell, normalIn_Py[1], M);
  children = faceTree[face]->children;
  for (E_Int i = 0; i < 4; i++) LFT[i] = children[i]; 
  T6_get_ordered_data(M, NODES[0], reorient, LFT, local);
  assert(local[0] == NODES[0]);
  assert(local[1] == NODES[3]);
  assert(local[3] == NODES[8]);
  NODES[4] = local[2];
  NODES[10] = local[4];
  NODES[11] = local[5];

  // RGT (Out)
  face = pf[2];
  reorient = get_reorient(face, cell, normalIn_Py[2], M);
  children = faceTree[face]->children;
  for (E_Int i = 0; i < 4; i++) RGT[i] = children[i]; 
  T6_get_ordered_data(M, NODES[1], reorient, RGT, local);
  assert(local[0] == NODES[1]);
  assert(local[1] == NODES[2]);
  assert(local[2] == NODES[4]);
  assert(local[3] == NODES[6]);
  NODES[12] = local[4];
  NODES[13] = local[5];

  // FRO (In)
  face = pf[3];
  reorient = get_reorient(face, cell, normalIn_Py[3], M);
  children = faceTree[face]->children;
  for (E_Int i = 0; i < 4; i++) FRO[i] = children[i]; 
  T6_get_ordered_data(M, NODES[1], reorient, FRO, local);
  assert(local[0] == NODES[1]);
  assert(local[1] == NODES[0]);
  assert(local[2] == NODES[4]);
  assert(local[3] == NODES[5]);
  assert(local[4] == NODES[11]);
  assert(local[5] == NODES[13]);

  // BCK (Out)
  face = pf[4];
  reorient = get_reorient(face, cell, normalIn_Py[4], M);
  children = faceTree[face]->children;
  for (E_Int i = 0; i < 4; i++) BCK[i] = children[i]; 
  T6_get_ordered_data(M, NODES[2], reorient, BCK, local);
  assert(local[0] == NODES[2]);
  assert(local[1] == NODES[3]);
  assert(local[2] == NODES[4]);
  assert(local[3] == NODES[7]);
  assert(local[4] == NODES[10]);
  assert(local[5] == NODES[12]);
  
  // Set internal faces in ngon
  E_Int *ptr = &M->indPG[M->nfaces];

  // nfaces (RGT of ncells)
  ptr[1] = ptr[0] + 3;
  M->ngon[*ptr  ] = NODES[5];
  M->ngon[*ptr+1] = NODES[9];
  M->ngon[*ptr+2] = NODES[11];
  ptr++;

  // nfaces+1 (BCK of ncells)
  ptr[1] = ptr[0] + 3;
  M->ngon[*ptr  ] = NODES[9];
  M->ngon[*ptr+1] = NODES[8];
  M->ngon[*ptr+2] = NODES[11];
  ptr++;

  // nfaces+2 (LFT of ncells+1)
  ptr[1] = ptr[0] + 3;
  M->ngon[*ptr  ] = NODES[5];
  M->ngon[*ptr+1] = NODES[9];
  M->ngon[*ptr+2] = NODES[13];
  ptr++;

  // nfaces+3 (BCK of ncells+1)
  ptr[1] = ptr[0] + 3;
  M->ngon[*ptr  ] = NODES[6];
  M->ngon[*ptr+1] = NODES[9];
  M->ngon[*ptr+2] = NODES[13];
  ptr++;

  // nfaces+4 (LFT of ncells+2)
  ptr[1] = ptr[0] + 3;
  M->ngon[*ptr  ] = NODES[9];
  M->ngon[*ptr+1] = NODES[7];
  M->ngon[*ptr+2] = NODES[12];
  ptr++;

  // nfaces+5 (FRO of ncells+2)
  ptr[1] = ptr[0] + 3;
  M->ngon[*ptr  ] = NODES[6];
  M->ngon[*ptr+1] = NODES[9];
  M->ngon[*ptr+2] = NODES[12];
  ptr++;

  // nfaces+6 (RGT of ncells+3)
  ptr[1] = ptr[0] + 3;
  M->ngon[*ptr  ] = NODES[9];
  M->ngon[*ptr+1] = NODES[7];
  M->ngon[*ptr+2] = NODES[10];
  ptr++;

  // nfaces+7 (FRO of ncells+3)
  ptr[1] = ptr[0] + 3;
  M->ngon[*ptr  ] = NODES[9];
  M->ngon[*ptr+1] = NODES[8];
  M->ngon[*ptr+2] = NODES[10];
  ptr++;

  // nfaces+8 (BOT of ncells+4)
  ptr[1] = ptr[0] + 4;
  M->ngon[*ptr  ] = NODES[11];
  M->ngon[*ptr+1] = NODES[13];
  M->ngon[*ptr+2] = NODES[12];
  M->ngon[*ptr+3] = NODES[10];
  ptr++;

  // nfaces+9 (LFT of ncells+5)
  ptr[1] = ptr[0] + 3;
  M->ngon[*ptr  ] = NODES[11];
  M->ngon[*ptr+1] = NODES[13];
  M->ngon[*ptr+2] = NODES[9];
  ptr++;

  // nfaces+10 (RGT of ncells+5)
  ptr[1] = ptr[0] + 3;
  M->ngon[*ptr  ] = NODES[10];
  M->ngon[*ptr+1] = NODES[12];
  M->ngon[*ptr+2] = NODES[9];
  ptr++;

  // nfaces+11 (FRO of ncells+5)
  ptr[1] = ptr[0] + 3;
  M->ngon[*ptr  ] = NODES[10];
  M->ngon[*ptr+1] = NODES[11];
  M->ngon[*ptr+2] = NODES[9];
  ptr++;

  // nfaces+12 (BCK of ncells+5)
  ptr[1] = ptr[0] + 3;
  M->ngon[*ptr  ] = NODES[12];
  M->ngon[*ptr+1] = NODES[13];
  M->ngon[*ptr+2] = NODES[9];
  ptr++;

  // Assemble children

  // First child replaces cell
  ptr = &M->indPH[cell];
  M->nface[*ptr  ] = BOT[0];
  M->nface[*ptr+1] = LFT[0];
  M->nface[*ptr+2] = M->nfaces;
  M->nface[*ptr+3] = FRO[1];
  M->nface[*ptr+4] = M->nfaces+1;

  ptr = &M->indPH[M->ncells];

  // ncells
  ptr[1] = ptr[0] + 5;
  M->nface[*ptr  ] = BOT[1];
  M->nface[*ptr+1] = M->nfaces+2;
  M->nface[*ptr+2] = RGT[0];
  M->nface[*ptr+3] = FRO[0];
  M->nface[*ptr+4] = M->nfaces+3;
  ptr++;

  // ncells+1
  ptr[1] = ptr[0] + 5;
  M->nface[*ptr  ] = BOT[2];
  M->nface[*ptr+1] = M->nfaces+4;
  M->nface[*ptr+2] = RGT[1];
  M->nface[*ptr+3] = M->nfaces+5;
  M->nface[*ptr+4] = BCK[0];
  ptr++;

  // ncells+2
  ptr[1] = ptr[0] + 5;
  M->nface[*ptr  ] = BOT[3];
  M->nface[*ptr+1] = LFT[1];
  M->nface[*ptr+2] = M->nfaces+6;
  M->nface[*ptr+3] = M->nfaces+7;
  M->nface[*ptr+4] = BCK[1];
  ptr++;

  // ncells+3
  ptr[1] = ptr[0] + 5;
  M->nface[*ptr  ] = M->nfaces+8;
  M->nface[*ptr+1] = LFT[2];
  M->nface[*ptr+2] = RGT[2];
  M->nface[*ptr+3] = FRO[2];
  M->nface[*ptr+4] = BCK[2];
  ptr++;

  // ncells+4
  ptr[1] = ptr[0] + 5;
  M->nface[*ptr  ] = M->nfaces+8;
  M->nface[*ptr+1] = M->nfaces+9;
  M->nface[*ptr+2] = M->nfaces+10;
  M->nface[*ptr+3] = M->nfaces+11;
  M->nface[*ptr+4] = M->nfaces+12;
  ptr++;

  // ncells+5
  ptr[1] = ptr[0] + 4;
  M->nface[*ptr  ] = M->nfaces;
  M->nface[*ptr+1] = FRO[3];
  M->nface[*ptr+2] = M->nfaces+9;
  M->nface[*ptr+3] = M->nfaces+2;
  ptr++;

  // ncells+6
  ptr[1] = ptr[0] + 4;
  M->nface[*ptr  ] = M->nfaces+3;
  M->nface[*ptr+1] = RGT[3];
  M->nface[*ptr+2] = M->nfaces+12;
  M->nface[*ptr+3] = M->nfaces+5;
  ptr++;

  // ncells+7
  ptr[1] = ptr[0] + 4;
  M->nface[*ptr  ] = M->nfaces+4;
  M->nface[*ptr+1] = M->nfaces+6;
  M->nface[*ptr+2] = BCK[3];
  M->nface[*ptr+3] = M->nfaces+10;
  ptr++;

  // ncells+8
  ptr[1] = ptr[0] + 4;
  M->nface[*ptr  ] = M->nfaces+11;
  M->nface[*ptr+1] = LFT[3];
  M->nface[*ptr+2] = M->nfaces+1;
  M->nface[*ptr+3] = M->nfaces+7;

  // Set new cells in tree
  E_Int next_level = M->cellTree[cell]->level + 1;

  // Parent
  Element *Elem = M->cellTree[cell];
  Elem->children = (E_Int *)XMALLOC(10 * sizeof(E_Int));
  Elem->children[0] = cell;
  for (E_Int i = 1; i < 10; i++) Elem->children[i] = M->ncells+i-1;
  Elem->nchildren = 10;
  Elem->parent = cell;
  Elem->position = 0;
  Elem->type = PYRA;
  Elem->level = next_level;
  Elem->state = 1;

  // Pyra children
  for (E_Int i = 0; i < 5; i++) {
    Element *Elem = M->cellTree[M->ncells+i];
    Elem->children = NULL;
    Elem->nchildren = 0;
    Elem->parent = cell;
    Elem->position = i+1;
    Elem->type = PYRA;
    Elem->level = next_level;
    Elem->state = 0;
  }

  // Tetra children
  for (E_Int i = 5; i < 9; i++) {
    Element *Elem = M->cellTree[M->ncells+i];
    Elem->children = NULL;
    Elem->nchildren = 0;
    Elem->parent = cell;
    Elem->position = i+1;
    Elem->type = TETRA;
    Elem->level = next_level;
    Elem->state = 0;
  }

  //set_parent_elem(M->cellTree, cell, 10, PYRA, M->ncells);
  //for (E_Int i = 0; i < 5; i++)
  //  set_child_elem(M->cellTree, cell, i+1, PYRA, next_level, M->ncells);
  //for (E_Int i = 5; i < 9; i++)
  //  set_child_elem(M->cellTree, cell, i+1, TETRA, next_level, M->ncells);

  // Set external faces owns and neis
  update_external_own_nei(cell, M);

  // tetras: ncells+6, ncells+7, ncells+8, ncells+9

  // Set owns and neis of internal faces
  M->owner[M->nfaces]     = cell;   
  M->neigh[M->nfaces]     = M->ncells+5;       
  M->owner[M->nfaces+1]   = cell;
  M->neigh[M->nfaces+1]   = M->ncells+8;

  M->owner[M->nfaces+2]   = M->ncells+5; 
  M->neigh[M->nfaces+2]   = M->ncells; 
  M->owner[M->nfaces+3]   = M->ncells;
  M->neigh[M->nfaces+3]   = M->ncells+6;
  
  M->owner[M->nfaces+4]   = M->ncells+7; 
  M->neigh[M->nfaces+4]   = M->ncells+1; 
  M->owner[M->nfaces+5]   = M->ncells+6;
  M->neigh[M->nfaces+5]   = M->ncells+1;
  
  M->owner[M->nfaces+6]   = M->ncells+2;
  M->neigh[M->nfaces+6]   = M->ncells+7;
  M->owner[M->nfaces+7]   = M->ncells+8;
  M->neigh[M->nfaces+7]   = M->ncells+2;

  M->owner[M->nfaces+8]   = M->ncells+4;
  M->neigh[M->nfaces+8]   = M->ncells+3;

  M->owner[M->nfaces+9]   = M->ncells+5;
  M->neigh[M->nfaces+9]   = M->ncells+4;
  M->owner[M->nfaces+10]  = M->ncells+4;
  M->neigh[M->nfaces+10]  = M->ncells+7;
  M->owner[M->nfaces+11]  = M->ncells+8;
  M->neigh[M->nfaces+11]  = M->ncells+4;
  M->owner[M->nfaces+12]  = M->ncells+4;
  M->neigh[M->nfaces+12]  = M->ncells+6;

  check_canon_pyra(cell, M);

  for (E_Int i = 0; i < 5; i++)
    check_canon_pyra(M->ncells+i, M);

  for (E_Int i = 5; i < 9; i++)
    check_canon_tetra(M->ncells+i, M);

  M->ncells += 9;
  M->nfaces += 13; 
  */
}

void reorder_pyra(E_Int i, E_Int nf, E_Int *pf, AMesh *M)
{
  E_Int common[4], map[4];
  E_Int bot = pf[0];
  E_Int *pn = get_facets(bot, M->ngon, M->indPG);
  
  for (E_Int i = 0; i < 4; i++) map[i] = pn[i];
  E_Int reorient = get_reorient(bot, i, normalIn_Py[0], M);
  if (reorient) std::swap(map[1], map[3]);

  E_Int lft, rgt, fro, bck;
  lft = rgt = fro = bck = -1;

  for (E_Int j = 1; j < 5; j++) {
    for (E_Int k = 0; k < 4; k++) common[k] = 0;

    E_Int face = pf[j];
    E_Int *pn = get_facets(face, M->ngon, M->indPG);

    for (E_Int k = 0; k < 3; k++) {
      E_Int point = pn[k];

      for (E_Int l = 0; l < 4; l++) {
        if (map[l] == point) {
          common[l] = 1;
          break;
        }
      }
    }

    if      (common[0] && common[3]) lft = face;
    else if (common[1] && common[2]) rgt = face;
    else if (common[1] && common[0]) fro = face;
    else                             bck = face;
  }

  assert(lft != -1);
  assert(rgt != -1);
  assert(fro != -1);
  assert(bck != -1);

  pf[1] = lft;
  pf[2] = rgt;
  pf[3] = fro;
  pf[4] = bck;
}

E_Int check_canon_pyra(E_Int cell, AMesh *M)
{
  E_Int NODES[5] = {-1, -1, -1, -1, -1};

  E_Int *pf = get_facets(cell, M->nface, M->indPH);

  E_Int face, i0, reorient, *pn, local[4];

  // BOT (in)
  face = pf[0];
  pn = get_facets(face, M->ngon, M->indPG);
  i0 = 0;
  reorient = get_reorient(face, cell, normalIn_Py[0], M);
  reorder_quad(local, pn, reorient, i0);
  for (E_Int i = 0; i < 4; i++) NODES[i] = local[i];

  // LFT (in)
  face = pf[1];
  pn = get_facets(face, M->ngon, M->indPG);
  i0 = Get_pos(NODES[0], pn, 3);
  reorient = get_reorient(face, cell, normalIn_Py[1], M);
  reorder_tri(local, pn, reorient, i0);
  assert(local[0] == NODES[0]);
  assert(local[1] == NODES[3]);
  NODES[4] = local[2];

  // RGT (out)
  face = pf[2];
  pn = get_facets(face, M->ngon, M->indPG);
  i0 = Get_pos(NODES[1], pn, 3);
  reorient = get_reorient(face, cell, normalIn_Py[2], M);
  reorder_tri(local, pn, reorient, i0);
  assert(local[0] == NODES[1]);
  assert(local[1] == NODES[2]);
  assert(local[2] == NODES[4]);

  // FRO (in)
  face = pf[3];
  pn = get_facets(face, M->ngon, M->indPG);
  i0 = Get_pos(NODES[1], pn, 3);
  reorient = get_reorient(face, cell, normalIn_Py[3], M);
  reorder_tri(local, pn, reorient, i0);
  assert(local[0] == NODES[1]);
  assert(local[1] == NODES[0]);
  assert(local[2] == NODES[4]);

  // BCK (out)
  face = pf[4];
  pn = get_facets(face, M->ngon, M->indPG);
  i0 = Get_pos(NODES[2], pn, 3);
  reorient = get_reorient(face, cell, normalIn_Py[4], M);
  reorder_tri(local, pn, reorient, i0);
  assert(local[0] == NODES[2]);
  assert(local[1] == NODES[3]);
  assert(local[2] == NODES[4]);

  return 1;
}

void make_ref_data_pyra(E_Int cell, AMesh *M, E_Float *pM, const pDirs &Dirs)
{}

void make_pdirs_pyra(E_Int cell, AMesh *M, pDirs &Dirs)
{}