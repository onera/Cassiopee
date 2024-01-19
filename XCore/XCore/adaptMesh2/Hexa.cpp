#include "Proto.h"

void refine_hexa(E_Int cell, AMesh *M)
{
  assert(M->cellTree->state(cell) != GONE);

  E_Int NODES[27], FACES[24];
  for (E_Int i = 0; i < 27; i++) NODES[i] = -1;
  E_Int *BOT = FACES;
  E_Int *TOP = FACES + 4;
  E_Int *LFT = FACES + 8;
  E_Int *RGT = FACES + 12;
  E_Int *FRO = FACES + 16;
  E_Int *BCK = FACES + 20;
  E_Int *pf = get_facets(cell, M->nface, M->indPH);

  E_Int face, reorient, local[9];
  Tree *faceTree = M->faceTree;
  Children *children;

  // BOT
  face = pf[0];
  reorient = get_reorient(face, cell, normalIn_H[0], M);
  children = faceTree->children(face);
  for (E_Int i = 0; i < 4; i++) BOT[i] = children->pc[i];
  NODES[0] = get_facets(face, M->ngon, M->indPG)[0];
  Q9_get_ordered_data(M, NODES[0], reorient, BOT, local);
  assert(local[0] == NODES[0]);
  for (E_Int i = 1; i < 4; i++) NODES[i] = local[i];
  NODES[8] = local[4];
  NODES[9] = local[5];
  NODES[10] = local[6];
  NODES[11] = local[7];
  NODES[12] = local[8];

  // LFT
  face = pf[2];
  reorient = get_reorient(face, cell, normalIn_H[2], M);
  children = faceTree->children(face);
  for (E_Int i = 0; i < 4; i++) LFT[i] = children->pc[i];
  Q9_get_ordered_data(M, NODES[0], reorient, LFT, local);
  assert(local[0] == NODES[0]);
  assert(local[1] == NODES[3]);
  assert(local[4] == NODES[11]);
  NODES[7] = local[2];
  NODES[4] = local[3];
  NODES[13] = local[5];
  NODES[14] = local[6];
  NODES[15] = local[7];
  NODES[16] = local[8];

  // RGT
  face = pf[3];
  reorient = get_reorient(face, cell, normalIn_H[3], M);
  children = faceTree->children(face);
  for (E_Int i = 0; i < 4; i++) RGT[i] = children->pc[i];
  Q9_get_ordered_data(M, NODES[1], reorient, RGT, local);
  assert(local[0] == NODES[1]);
  assert(local[1] == NODES[2]);
  assert(local[4] == NODES[9]);
  NODES[6] = local[2];
  NODES[5] = local[3];
  NODES[17] = local[5];
  NODES[18] = local[6];
  NODES[19] = local[7];
  NODES[20] = local[8];

  // FRO
  face = pf[4];
  reorient = get_reorient(face, cell, normalIn_H[4], M);
  children = faceTree->children(face);
  for (E_Int i = 0; i < 4; i++) FRO[i] = children->pc[i];
  Q9_get_ordered_data(M, NODES[1], reorient, FRO, local);
  assert(local[0] == NODES[1]);
  assert(local[1] == NODES[0]);
  assert(local[2] == NODES[4]);
  assert(local[3] == NODES[5]);
  assert(local[4] == NODES[8]);
  assert(local[5] == NODES[15]);
  assert(local[7] == NODES[19]);
  NODES[21] = local[6];
  NODES[22] = local[8];

  // BCK
  face = pf[5];
  reorient = get_reorient(face, cell, normalIn_H[5], M);
  children = faceTree->children(face);
  for (E_Int i = 0; i < 4; i++) BCK[i] = children->pc[i];
  Q9_get_ordered_data(M, NODES[2], reorient, BCK, local);
  assert(local[0] == NODES[2]);
  assert(local[1] == NODES[3]);
  assert(local[2] == NODES[7]);
  assert(local[3] == NODES[6]);
  assert(local[4] == NODES[10]);
  assert(local[5] == NODES[13]);
  assert(local[7] == NODES[17]);
  NODES[23] = local[6];
  NODES[24] = local[8];

  // TOP
  face = pf[1];
  reorient = get_reorient(face, cell, normalIn_H[1], M);
  children = faceTree->children(face);
  for (E_Int i = 0; i < 4; i++) TOP[i] = children->pc[i];
  Q9_get_ordered_data(M, NODES[4], reorient, TOP, local);
  assert(local[0] == NODES[4]);
  assert(local[1] == NODES[5]);
  assert(local[2] == NODES[6]);
  assert(local[3] == NODES[7]);
  assert(local[4] == NODES[21]);
  assert(local[5] == NODES[18]);
  assert(local[6] == NODES[23]);
  assert(local[7] == NODES[14]);
  NODES[25] = local[8];

  // Add cell centroid
  // Note(Imad): supposed to be already computed
  NODES[26] = M->npoints;
  M->x[M->npoints] = M->cx[cell];
  M->y[M->npoints] = M->cy[cell];
  M->z[M->npoints] = M->cz[cell];
  M->npoints++;

  // Set internal faces in ngon
  E_Int *ptr = &M->indPG[M->nfaces];

  // NCELLS TOP && RGT && BCK

  // nfaces
  ptr[1] = ptr[0] + 4;
  M->ngon[*ptr  ] = NODES[15]; M->ngon[*ptr+1] = NODES[22];
  M->ngon[*ptr+2] = NODES[26]; M->ngon[*ptr+3] = NODES[16];
  ptr++;

  // nfaces+1
  ptr[1] = ptr[0] + 4;
  M->ngon[*ptr  ] = NODES[8];  M->ngon[*ptr+1] = NODES[12];
  M->ngon[*ptr+2] = NODES[26]; M->ngon[*ptr+3] = NODES[22];
  ptr++;

  // nfaces+2
  ptr[1] = ptr[0] + 4;
  M->ngon[*ptr  ] = NODES[12]; M->ngon[*ptr+1] = NODES[11];
  M->ngon[*ptr+2] = NODES[16]; M->ngon[*ptr+3] = NODES[26];
  ptr++;

  // NCELLS+1 TOP && BCK

  // nfaces+3
  ptr[1] = ptr[0] + 4;
  M->ngon[*ptr  ] = NODES[22]; M->ngon[*ptr+1] = NODES[19];
  M->ngon[*ptr+2] = NODES[20]; M->ngon[*ptr+3] = NODES[26];
  ptr++;

  // nfaces+4
  ptr[1] = ptr[0] + 4;
  M->ngon[*ptr  ] = NODES[9]; M->ngon[*ptr+1] = NODES[12];
  M->ngon[*ptr+2] = NODES[26]; M->ngon[*ptr+3] = NODES[20];
  ptr++;

  // NCELLS+2 TOP && LFT

  // nfaces+5
  ptr[1] = ptr[0] + 4;
  M->ngon[*ptr  ] = NODES[26]; M->ngon[*ptr+1] = NODES[20];
  M->ngon[*ptr+2] = NODES[17]; M->ngon[*ptr+3] = NODES[24];
  ptr++;

  // nfaces+6
  ptr[1] = ptr[0] + 4;
  M->ngon[*ptr  ] = NODES[12]; M->ngon[*ptr+1] = NODES[10];
  M->ngon[*ptr+2] = NODES[24]; M->ngon[*ptr+3] = NODES[26];
  ptr++;

  // NCELLS+3 TOP

  // nfaces+7
  ptr[1] = ptr[0] + 4;
  M->ngon[*ptr  ] = NODES[16];  M->ngon[*ptr+1] = NODES[26];
  M->ngon[*ptr+2] = NODES[24]; M->ngon[*ptr+3] = NODES[13];
  ptr++;

  // NCELLS+4 RGT && BCK

  // nfaces+8
  ptr[1] = ptr[0] + 4;
  M->ngon[*ptr  ] = NODES[22]; M->ngon[*ptr+1] = NODES[26];
  M->ngon[*ptr+2] = NODES[25]; M->ngon[*ptr+3] = NODES[21];
  ptr++;

  // nfaces+9
  ptr[1] = ptr[0] + 4;
  M->ngon[*ptr  ] = NODES[26]; M->ngon[*ptr+1] = NODES[16];
  M->ngon[*ptr+2] = NODES[14]; M->ngon[*ptr+3] = NODES[25];
  ptr++;

  // NCELLS+5 BCK

  // nfaces+10
  ptr[1] = ptr[0] + 4;
  M->ngon[*ptr  ] = NODES[20]; M->ngon[*ptr+1] = NODES[26];
  M->ngon[*ptr+2] = NODES[25]; M->ngon[*ptr+3] = NODES[18];
  ptr++;

  // NCELLS+6 LFT

  // nfaces+11
  ptr[1] = ptr[0] + 4;
  M->ngon[*ptr  ] = NODES[26]; M->ngon[*ptr+1] = NODES[24];
  M->ngon[*ptr+2] = NODES[23]; M->ngon[*ptr+3] = NODES[25];

  // Assemble children

  // First child replaces cell
  ptr = &M->indPH[cell];
  M->nface[*ptr  ] = BOT[0];       M->nface[*ptr+1] = M->nfaces;
  M->nface[*ptr+2] = LFT[0];       M->nface[*ptr+3] = M->nfaces+1;
  M->nface[*ptr+4] = FRO[1];       M->nface[*ptr+5] = M->nfaces+2;

  ptr = &M->indPH[M->ncells];

  // ncells
  ptr[1] = ptr[0] + 6;
  M->nface[*ptr  ] = BOT[1];       M->nface[*ptr+1] = M->nfaces+3;
  M->nface[*ptr+2] = M->nfaces+1;  M->nface[*ptr+3] = RGT[0];
  M->nface[*ptr+4] = FRO[0];       M->nface[*ptr+5] = M->nfaces+4;
  ptr++;

  // ncells+1
  ptr[1] = ptr[0] + 6;
  M->nface[*ptr  ] = BOT[2];       M->nface[*ptr+1] = M->nfaces+5;
  M->nface[*ptr+2] = M->nfaces+6;  M->nface[*ptr+3] = RGT[1];
  M->nface[*ptr+4] = M->nfaces+4;  M->nface[*ptr+5] = BCK[0];
  ptr++;

  // ncells+2
  ptr[1] = ptr[0] + 6;
  M->nface[*ptr  ] = BOT[3];       M->nface[*ptr+1] = M->nfaces+7;
  M->nface[*ptr+2] = LFT[1];       M->nface[*ptr+3] = M->nfaces+6;
  M->nface[*ptr+4] = M->nfaces+2;  M->nface[*ptr+5] = BCK[1];
  ptr++;

  /*********/

  // ncells+3
  ptr[1] = ptr[0] + 6;
  M->nface[*ptr  ] = M->nfaces;    M->nface[*ptr+1] = TOP[0];
  M->nface[*ptr+2] = LFT[3];       M->nface[*ptr+3] = M->nfaces+8;
  M->nface[*ptr+4] = FRO[2];       M->nface[*ptr+5] = M->nfaces+9;
  ptr++;

  // ncells+4
  ptr[1] = ptr[0] + 6;
  M->nface[*ptr  ] = M->nfaces+3;  M->nface[*ptr+1] = TOP[1];
  M->nface[*ptr+2] = M->nfaces+8;  M->nface[*ptr+3] = RGT[3];
  M->nface[*ptr+4] = FRO[3];       M->nface[*ptr+5] = M->nfaces+10;
  ptr++;

  // ncells+5
  ptr[1] = ptr[0] + 6;
  M->nface[*ptr  ] = M->nfaces+5; M->nface[*ptr+1] = TOP[2];
  M->nface[*ptr+2] = M->nfaces+11; M->nface[*ptr+3] = RGT[2];
  M->nface[*ptr+4] = M->nfaces+10; M->nface[*ptr+5] = BCK[3];
  ptr++;

  // ncells+6
  ptr[1] = ptr[0] + 6;
  M->nface[*ptr  ] = M->nfaces+7;  M->nface[*ptr+1] = TOP[3];
  M->nface[*ptr+2] = LFT[2];  M->nface[*ptr+3] = M->nfaces+11;
  M->nface[*ptr+4] = M->nfaces+9; M->nface[*ptr+5] = BCK[2];

  // Set new cells in tree

  E_Int next_level = M->cellTree->level(cell) + 1;

  // Parent
  M->cellTree->set_parent_elem(cell, 8, M->ncells);

  // Hexa children
  for (E_Int i = 1; i < 8; i++) {
    M->cellTree->set_child_elem(i, cell, HEXA, next_level, M->ncells);
  }

  // Set external faces owns and neis
  update_external_pe_after_ref(cell, M);

  // Set internal faces in faceTree
  for (E_Int i = 0; i < 12; i++)
    M->faceTree->set_new_elem(M->nfaces+i, QUAD, next_level);

  // Set owns and neis of internal faces
  M->owner[M->nfaces]    = cell;        M->owner[M->nfaces+1]  = cell;
  M->neigh[M->nfaces]    = M->ncells+3; M->neigh[M->nfaces+1]  = M->ncells  ;

  M->owner[M->nfaces+2]  = cell;        M->owner[M->nfaces+3]  = M->ncells  ;
  M->neigh[M->nfaces+2]  = M->ncells+2; M->neigh[M->nfaces+3]  = M->ncells+4;

  M->owner[M->nfaces+4]  = M->ncells  ; M->owner[M->nfaces+5]  = M->ncells+1;
  M->neigh[M->nfaces+4]  = M->ncells+1; M->neigh[M->nfaces+5]  = M->ncells+5;

  M->owner[M->nfaces+6]  = M->ncells+2; M->owner[M->nfaces+7]  = M->ncells+2;
  M->neigh[M->nfaces+6]  = M->ncells+1; M->neigh[M->nfaces+7]  = M->ncells+6;

  M->owner[M->nfaces+8]  = M->ncells+3; M->owner[M->nfaces+9]  = M->ncells+3;
  M->neigh[M->nfaces+8]  = M->ncells+4; M->neigh[M->nfaces+9]  = M->ncells+6;

  M->owner[M->nfaces+10] = M->ncells+4; M->owner[M->nfaces+11] = M->ncells+6;
  M->neigh[M->nfaces+10] = M->ncells+5; M->neigh[M->nfaces+11] = M->ncells+5;

  check_canon_hexa(cell, M);
  for (E_Int i = 0; i < 7; i++)
    check_canon_hexa(M->ncells+i, M);

  M->ncells += 7;
  M->nfaces += 12;
}
