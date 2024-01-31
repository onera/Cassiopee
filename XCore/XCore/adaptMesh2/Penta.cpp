#include "Proto.h"

void Pe18_refine(E_Int cell, AMesh *M)
{
  /*
  E_Int NODES[18], FACES[20];
  for (E_Int i = 0; i < 18; i++) NODES[i] = -1;
  E_Int *BOT = FACES;
  E_Int *TOP = FACES + 4;
  E_Int *LFT = FACES + 8;
  E_Int *RGT = FACES + 12;
  E_Int *BCK = FACES + 16;
  E_Int *pf = get_facets(cell, M->nface, M->indPH);

  E_Int face, i0, reorient, *children, nchildren, local[9];
  Element **faceTree = M->faceTree;

  // BOT
  face = pf[0];
  reorient = get_reorient(face, cell, normalIn_Pe[0], M);
  children = faceTree[face]->children;
  for (E_Int i = 0; i < 4; i++) BOT[i] = children[i];
  NODES[0] = get_facets(face, M->ngon, M->indPG)[0];
  T6_get_ordered_data(M, NODES[0], reorient, BOT, local);
  assert(NODES[0] == local[0]);
  for (E_Int i = 1; i < 3; i++) NODES[i] = local[i];
  NODES[6] = local[3];
  NODES[7] = local[4];
  NODES[8] = local[5];

  // LFT
  face = pf[2];
  reorient = get_reorient(face, cell, normalIn_Pe[2], M);
  children = faceTree[face]->children;
  for (E_Int i = 0; i < 4; i++) LFT[i] = children[i]; 
  Q9_get_ordered_data(M, NODES[0], reorient, LFT, local);
  assert(local[0] == NODES[0]);
  assert(local[1] == NODES[2]);
  assert(local[4] == NODES[8]);
  NODES[5] = local[2];
  NODES[3] = local[3];
  NODES[9] = local[5];
  NODES[10] = local[6];
  NODES[11] = local[7];
  NODES[12] = local[8];

  // RGT
  face = pf[3];
  reorient = get_reorient(face, cell, normalIn_Pe[3], M);
  children = faceTree[face]->children;
  for (E_Int i = 0; i < 4; i++) RGT[i] = children[i]; 
  Q9_get_ordered_data(M, NODES[0], reorient, RGT, local);
  assert(local[0] == NODES[0]);
  assert(local[1] == NODES[1]);
  assert(local[3] == NODES[3]);
  assert(local[4] == NODES[6]);
  assert(local[7] == NODES[11]);
  NODES[4] = local[2];
  NODES[13] = local[5];
  NODES[14] = local[6];
  NODES[15] = local[8];

  // BCK (In)
  face = pf[4];
  reorient = get_reorient(face, cell, normalIn_Pe[4], M);
  children = faceTree[face]->children;
  for (E_Int i = 0; i < 4; i++) BCK[i] = children[i]; 
  Q9_get_ordered_data(M, NODES[2], reorient, BCK, local);
  assert(local[0] == NODES[2]);
  assert(local[1] == NODES[1]);
  assert(local[2] == NODES[4]);
  assert(local[3] == NODES[5]);
  assert(local[4] == NODES[7]);
  assert(local[5] == NODES[13]);
  assert(local[7] == NODES[9]);
  NODES[16] = local[6];
  NODES[17] = local[8];

  // TOP
  face = pf[1];
  reorient = get_reorient(face, cell, normalIn_Pe[1], M);
  children = faceTree[face]->children;
  for (E_Int i = 0; i < 4; i++) TOP[i] = children[i]; 
  T6_get_ordered_data(M, NODES[3], reorient, TOP, local);
  assert(local[0] == NODES[3]);
  assert(local[1] == NODES[4]);
  assert(local[2] == NODES[5]);
  assert(local[3] == NODES[14]);
  assert(local[4] == NODES[16]);
  assert(local[5] == NODES[10]);

  // Set internal faces in ngon
  E_Int *ptr = &M->indPG[M->nfaces];

  // nfaces (TOP of ncells)
  ptr[1] = ptr[0] + 3;
  M->ngon[*ptr  ] = NODES[11];
  M->ngon[*ptr+1] = NODES[15];
  M->ngon[*ptr+2] = NODES[12];
  ptr++;

  // nfaces+1 (BCK of ncells)
  ptr[1] = ptr[0] + 4;
  M->ngon[*ptr  ] = NODES[8];
  M->ngon[*ptr+1] = NODES[6];
  M->ngon[*ptr+2] = NODES[15];
  M->ngon[*ptr+3] = NODES[12];
  ptr++;

  // nfaces+2 (TOP of ncells+1)
  ptr[1] = ptr[0] + 3;
  M->ngon[*ptr  ] = NODES[15];
  M->ngon[*ptr+1] = NODES[13];
  M->ngon[*ptr+2] = NODES[17];
  ptr++;

  // nfaces+3 (LFT of ncells+1)
  ptr[1] = ptr[0] + 4;
  M->ngon[*ptr  ] = NODES[6];
  M->ngon[*ptr+1] = NODES[7];
  M->ngon[*ptr+2] = NODES[17];
  M->ngon[*ptr+3] = NODES[15];
  ptr++;

  // nfaces+4 (TOP of ncells+2)
  ptr[1] = ptr[0] + 3;
  M->ngon[*ptr  ] = NODES[12];
  M->ngon[*ptr+1] = NODES[17];
  M->ngon[*ptr+2] = NODES[9];
  ptr++;

  // nfaces+5 (RGT of ncells+2)
  ptr[1] = ptr[0] + 4;
  M->ngon[*ptr  ] = NODES[8];
  M->ngon[*ptr+1] = NODES[7];
  M->ngon[*ptr+2] = NODES[17];
  M->ngon[*ptr+3] = NODES[12];
  ptr++;

  // nfaces+6 (BOT of ncells+3)
  ptr[1] = ptr[0] + 3;
  M->ngon[*ptr  ] = NODES[17];
  M->ngon[*ptr+1] = NODES[15];
  M->ngon[*ptr+2] = NODES[12];
  ptr++;

  // nfaces+7 (BCK of ncells+4)
  ptr[1] = ptr[0] + 4;
  M->ngon[*ptr  ] = NODES[12];
  M->ngon[*ptr+1] = NODES[15];
  M->ngon[*ptr+2] = NODES[14];
  M->ngon[*ptr+3] = NODES[10];
  ptr++;

  // nfaces+8 (LFT of ncells+5)
  ptr[1] = ptr[0] + 4;
  M->ngon[*ptr  ] = NODES[15];
  M->ngon[*ptr+1] = NODES[17];
  M->ngon[*ptr+2] = NODES[16];
  M->ngon[*ptr+3] = NODES[14];
  ptr++;

  // nfaces+9 (RGT of ncells+6)
  ptr[1] = ptr[0] + 4;
  M->ngon[*ptr  ] = NODES[12];
  M->ngon[*ptr+1] = NODES[17];
  M->ngon[*ptr+2] = NODES[16];
  M->ngon[*ptr+3] = NODES[10];

  // Assemble children

  // First child replaces cell
  ptr = &M->indPH[cell];
  M->nface[*ptr  ] = BOT[0];
  M->nface[*ptr+1] = M->nfaces;
  M->nface[*ptr+2] = LFT[0];
  M->nface[*ptr+3] = RGT[0];
  M->nface[*ptr+4] = M->nfaces+1;

  ptr = &M->indPH[M->ncells];

  // ncells
  ptr[1] = ptr[0] + 5;
  M->nface[*ptr  ] = BOT[1];
  M->nface[*ptr+1] = M->nfaces+2;
  M->nface[*ptr+2] = M->nfaces+3;
  M->nface[*ptr+3] = RGT[1];
  M->nface[*ptr+4] = BCK[1];
  ptr++;

  // ncells+1
  ptr[1] = ptr[0] + 5;
  M->nface[*ptr  ] = BOT[2];
  M->nface[*ptr+1] = M->nfaces+4;
  M->nface[*ptr+2] = LFT[1];
  M->nface[*ptr+3] = M->nfaces+5;
  M->nface[*ptr+4] = BCK[0];
  ptr++;

  // ncells+2
  ptr[1] = ptr[0] + 5;
  M->nface[*ptr  ] = M->nfaces+6;
  M->nface[*ptr+1] = BOT[3];
  M->nface[*ptr+2] = M->nfaces+5;
  M->nface[*ptr+3] = M->nfaces+3;
  M->nface[*ptr+4] = M->nfaces+1;
  ptr++;

  // ncells+3
  ptr[1] = ptr[0] + 5;
  M->nface[*ptr  ] = M->nfaces;
  M->nface[*ptr+1] = TOP[0];
  M->nface[*ptr+2] = LFT[3];
  M->nface[*ptr+3] = RGT[3];
  M->nface[*ptr+4] = M->nfaces+7;
  ptr++;

  // ncells+4
  ptr[1] = ptr[0] + 5;
  M->nface[*ptr  ] = M->nfaces+2;
  M->nface[*ptr+1] = TOP[1];
  M->nface[*ptr+2] = M->nfaces+8;
  M->nface[*ptr+3] = RGT[2];
  M->nface[*ptr+4] = BCK[2];
  ptr++;

  // ncells+5
  ptr[1] = ptr[0] + 5;
  M->nface[*ptr  ] = M->nfaces+4;
  M->nface[*ptr+1] = TOP[2];
  M->nface[*ptr+2] = LFT[2];
  M->nface[*ptr+3] = M->nfaces+9;
  M->nface[*ptr+4] = BCK[3];
  ptr++;

  // ncells+6
  ptr[1] = ptr[0] + 5;
  M->nface[*ptr  ] = TOP[3];
  M->nface[*ptr+1] = M->nfaces+6;
  M->nface[*ptr+2] = M->nfaces+9;
  M->nface[*ptr+3] = M->nfaces+8;
  M->nface[*ptr+4] = M->nfaces+7;

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
  Elem->type = PENTA;
  Elem->level = next_level;
  Elem->state = 1;

  // Penta children
  for (E_Int i = 0; i < 7; i++) {
    Element *Elem = M->cellTree[M->ncells+i];
    Elem->children = NULL;
    Elem->nchildren = 0;
    Elem->parent = cell;
    Elem->position = i+1;
    Elem->type = PENTA;
    Elem->level = next_level;
    Elem->state = 0;
  }

  //E_Int next_level = M->cellTree[cell]->level + 1;
  //set_parent_elem(M->cellTree, cell, 8, PENTA, M->ncells);
  //for (E_Int i = 0; i < 7; i++)
  //  set_child_elem(M->cellTree, cell, i, PENTA, next_level, M->ncells);

  // Set external faces owns and neis
  update_external_own_nei(cell, M);

  // Set internal faces in faceTree
  set_new_face_in_tree(faceTree, M->nfaces,   TRI,  next_level);
  set_new_face_in_tree(faceTree, M->nfaces+1, QUAD, next_level);
  set_new_face_in_tree(faceTree, M->nfaces+2, TRI,  next_level);
  set_new_face_in_tree(faceTree, M->nfaces+3, QUAD, next_level);
  set_new_face_in_tree(faceTree, M->nfaces+4, TRI,  next_level);
  set_new_face_in_tree(faceTree, M->nfaces+5, QUAD, next_level);
  set_new_face_in_tree(faceTree, M->nfaces+6, TRI,  next_level);
  set_new_face_in_tree(faceTree, M->nfaces+7, QUAD, next_level);
  set_new_face_in_tree(faceTree, M->nfaces+8, QUAD, next_level);
  set_new_face_in_tree(faceTree, M->nfaces+9, QUAD, next_level);

  // Set owns and neis of internal faces
  M->owner[M->nfaces]    = cell;   
  M->neigh[M->nfaces]    = M->ncells+3;       

  M->owner[M->nfaces+1]  = M->ncells+2;
  M->neigh[M->nfaces+1]  = cell;

  M->owner[M->nfaces+2]  = M->ncells; 
  M->neigh[M->nfaces+2]  = M->ncells+4; 

  M->owner[M->nfaces+3]  = M->ncells+2;
  M->neigh[M->nfaces+3]  = M->ncells;

  M->owner[M->nfaces+4]  = M->ncells+1; 
  M->neigh[M->nfaces+4]  = M->ncells+5; 

  M->owner[M->nfaces+5]  = M->ncells+1;
  M->neigh[M->nfaces+5]  = M->ncells+2;

  M->owner[M->nfaces+6]  = M->ncells+6;
  M->neigh[M->nfaces+6]  = M->ncells+2;

  M->owner[M->nfaces+7]  = M->ncells+6;
  M->neigh[M->nfaces+7]  = M->ncells+3;

  M->owner[M->nfaces+8]  = M->ncells+6;
  M->neigh[M->nfaces+8]  = M->ncells+4;

  M->owner[M->nfaces+9]  = M->ncells+5;
  M->neigh[M->nfaces+9]  = M->ncells+6;

  check_canon_penta(cell, M);
  for (E_Int i = 0; i < 7; i++)
    check_canon_penta(M->ncells+i, M);

  M->ncells += 7;
  M->nfaces += 10;
  */
}

void Pe18_unrefine(E_Int cell, AMesh *M)
{}
void Pe12_refine(E_Int cell, AMesh *M)
{}
void Pe12_unrefine(E_Int cell, AMesh *M)
{}

void reorder_penta(E_Int i, E_Int nf, E_Int *pf, AMesh *M)
{
  // First tri is bottom
  E_Int bot = -1;
  for (E_Int i = 0; i < nf; i++) {
    if (get_stride(pf[i], M->indPG) == 3) {
      bot = pf[i];
      Right_shift(pf, i, 5);
      break;
    }
  }
  assert(bot != -1);

  E_Int common[3], map[3];
  E_Int *pn = get_facets(bot, M->ngon, M->indPG);
  
  for (E_Int i = 0; i < 3; i++) map[i] = pn[i];
  E_Int reorient = get_reorient(bot, i, normalIn_Pe[0], M);
  if (reorient) std::swap(map[1], map[2]);

  E_Int top, lft, rgt, bck;
  top = lft = rgt = bck = -1;

  for (E_Int j = 1; j < 5; j++) {
    for (E_Int k = 0; k < 3; k++) common[k] = 0;

    E_Int face = pf[j];
    E_Int stride = get_stride(face, M->indPG);
    E_Int *pn = get_facets(face, M->ngon, M->indPG);

    for (E_Int k = 0; k < stride; k++) {
      E_Int point = pn[k];

      for (E_Int l = 0; l < 3; l++) {
        if (map[l] == point) {
          common[l] = 1;
          break;
        }
      }
    }

    if      (common[0] && common[2]) lft = face;
    else if (common[0] && common[1]) rgt = face;
    else if (common[1] && common[2]) bck = face;
    else                             top = face;
  }
  assert(top != -1);
  assert(lft != -1);
  assert(rgt != -1);
  assert(bck != -1);

  pf[1] = top;
  pf[2] = lft;
  pf[3] = rgt;
  pf[4] = bck;
}

E_Int check_canon_penta(E_Int cell, AMesh *M)
{
  E_Int NODES[6] = {-1, -1, -1, -1, -1, -1};

  E_Int *pf = get_facets(cell, M->nface, M->indPH);

  E_Int face, i0, reorient, *pn, local[4];

  // BOT
  face = pf[0];
  pn = get_facets(face, M->ngon, M->indPG);
  i0 = 0;
  reorient = get_reorient(face, cell, normalIn_Pe[0], M);
  reorder_tri(local, pn, reorient, i0);
  for (E_Int i = 0; i < 3; i++) NODES[i] = local[i];

  // LFT (in)
  face = pf[2];
  pn = get_facets(face, M->ngon, M->indPG);
  i0 = Get_pos(NODES[0], pn, 4);
  assert(i0 != -1);
  reorient = get_reorient(face, cell, normalIn_Pe[2], M);
  reorder_quad(local, pn, reorient, i0);
  assert(local[0] == NODES[0]);
  assert(local[1] == NODES[2]);
  NODES[5] = local[2];
  NODES[3] = local[3];

  // RGT (out)
  face = pf[3];
  pn = get_facets(face, M->ngon, M->indPG);
  i0 = Get_pos(NODES[0], pn, 4);
  reorient = get_reorient(face, cell, normalIn_Pe[3], M);
  reorder_quad(local, pn, reorient, i0);
  assert(local[0] == NODES[0]);
  assert(local[1] == NODES[1]);
  assert(local[3] == NODES[3]);
  NODES[4] = local[2];

  // BCK (in)
  face = pf[4];
  pn = get_facets(face, M->ngon, M->indPG);
  i0 = Get_pos(NODES[2], pn, 4);
  reorient = get_reorient(face, cell, normalIn_Pe[4], M);
  reorder_quad(local, pn, reorient, i0);
  assert(local[0] == NODES[2]);
  assert(local[1] == NODES[1]);
  assert(local[2] == NODES[4]);
  assert(local[3] == NODES[5]);

  // TOP
  face = pf[1];
  pn = get_facets(face, M->ngon, M->indPG);
  i0 = Get_pos(NODES[3], pn, 3);
  reorient = get_reorient(face, cell, normalIn_Pe[1], M);
  reorder_tri(local, pn, reorient, i0);
  assert(local[0] == NODES[3]);
  assert(local[1] == NODES[4]);
  assert(local[2] == NODES[5]);

  return 1;
}
