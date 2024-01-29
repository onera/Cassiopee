#include "Proto.h"

void H27_refine(E_Int cell, AMesh *M)
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

  assert(check_canon_hexa(cell, M));
  for (E_Int i = 0; i < 7; i++)
    assert(check_canon_hexa(M->ncells+i, M));

  M->ncells += 7;
  M->nfaces += 12;
}

void H27_unrefine(E_Int cell, AMesh *M)
{
  M->cellTree->state_[cell] = UNREFINED;
  M->cellTree->level_[cell] -= 1;

  Children *cur = M->cellTree->children(cell);

  assert(cur != NULL);

  for (E_Int i = 1; i < cur->n; i++) M->cellTree->state_[cur->pc[i]] = GONE;

  E_Int c0 = cur->pc[0];
  E_Int c1 = cur->pc[1];
  E_Int c2 = cur->pc[2];
  E_Int c3 = cur->pc[3];
  E_Int c4 = cur->pc[4];
  E_Int c5 = cur->pc[5];
  E_Int c6 = cur->pc[6];
  E_Int c7 = cur->pc[7];

  E_Int nf = -1;

  E_Int *f0 = get_cell(c0, nf, M->nface, M->indPH);
  E_Int *f1 = get_cell(c1, nf, M->nface, M->indPH);
  E_Int *f2 = get_cell(c2, nf, M->nface, M->indPH);
  E_Int *f3 = get_cell(c3, nf, M->nface, M->indPH);
  E_Int *f4 = get_cell(c4, nf, M->nface, M->indPH);
  E_Int *f5 = get_cell(c5, nf, M->nface, M->indPH);
  E_Int *f6 = get_cell(c6, nf, M->nface, M->indPH);
  E_Int *f7 = get_cell(c7, nf, M->nface, M->indPH);

  assert(f0[3] == f1[2]);
  assert(f1[5] == f2[4]);
  assert(f2[2] == f3[3]);
  assert(f3[4] == f0[5]);
  assert(f4[3] == f5[2]);
  assert(f5[5] == f6[4]);
  assert(f6[2] == f7[3]);
  assert(f7[4] == f4[5]);
  assert(f0[1] == f4[0]);
  assert(f1[1] == f5[0]);
  assert(f2[1] == f6[0]);
  assert(f3[1] == f7[0]);

  // Set owner/neigh/parent of internal faces to -1

  //printf("Gone faces: ");
  
  // rgt and top of c0
  M->faceTree->state_[f0[3]] = GONE;
  M->faceTree->state_[f0[1]] = GONE;
  M->owner[f0[3]] = -1;
  M->neigh[f0[3]] = -1;
  M->owner[f0[1]] = -1;
  M->neigh[f0[1]] = -1;

  // bck and top of c1
  M->faceTree->state_[f1[5]] = GONE;
  M->faceTree->state_[f1[1]] = GONE;
  M->owner[f1[5]] = -1;
  M->neigh[f1[5]] = -1;
  M->owner[f1[1]] = -1;
  M->neigh[f1[1]] = -1;
  
  // lft and top of c2
  M->faceTree->state_[f2[2]] = GONE;
  M->faceTree->state_[f2[1]] = GONE;
  M->owner[f2[2]] = -1;
  M->neigh[f2[2]] = -1;
  M->owner[f2[1]] = -1;
  M->neigh[f2[1]] = -1;

  // fro and top of c3
  M->faceTree->state_[f3[4]] = GONE;
  M->faceTree->state_[f3[1]] = GONE;
  M->owner[f3[4]] = -1;
  M->neigh[f3[4]] = -1;
  M->owner[f3[1]] = -1;
  M->neigh[f3[1]] = -1;

  // rgt of c4
  M->faceTree->state_[f4[3]] = GONE;
  M->owner[f4[3]] = -1;
  M->neigh[f4[3]] = -1;

  // bck of c5
  M->faceTree->state_[f5[5]] = GONE;
  M->owner[f5[5]] = -1;
  M->neigh[f5[5]] = -1;

  // lft of c6
  M->faceTree->state_[f6[2]] = GONE;
  M->owner[f6[2]] = -1;
  M->neigh[f6[2]] = -1;

  // fro of c7
  M->faceTree->state_[f7[4]] = GONE;
  M->owner[f7[4]] = -1;
  M->neigh[f7[4]] = -1;

  update_external_pe_after_unref(cell, M);

  // We know that these 6 faces are children of a prev refined face
  nf = -1;
  E_Int *pf = get_cell(cell, nf, M->nface, M->indPH);

  pf[0] = master_face(f0[0], M);
  pf[1] = master_face(f4[1], M);
  pf[2] = master_face(f0[2], M);
  pf[3] = master_face(f1[3], M);
  pf[4] = master_face(f0[4], M);
  pf[5] = master_face(f3[5], M);
}

void H18_refine(E_Int cell, AMesh *M)
{
  E_Int NODES[18], FACES[16];
  E_Int *BOT = FACES;
  E_Int *TOP = FACES + 4;
  E_Int *LFT = FACES + 8;
  E_Int *RGT = FACES + 10;
  E_Int *FRO = FACES + 12;
  E_Int *BCK = FACES + 14;
  E_Int *pf = get_facets(cell, M->nface, M->indPH);

  E_Int face, reorient, local[9];
  Children *children;
  Tree *faceTree = M->faceTree;

  // BOT
  face = pf[0];
  if (M->mode_2D)
    assert(check_face_aligned_with_vector(face, M->mode_2D, M));
  reorient = get_reorient(face, cell, normalIn_H[0], M);
  children = faceTree->children(face);
  assert(children->n == 4);
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
  assert(children->n == 2);
  for (E_Int i = 0; i < 2; i++) LFT[i] = children->pc[i];
  Q6_get_ordered_data(M, NODES[0], reorient, LFT, local);
  assert(local[0] == NODES[0]);
  assert(local[1] == NODES[3]);
  assert(local[4] == NODES[11]);
  NODES[7] = local[2];
  NODES[4] = local[3];
  NODES[13] = local[5];

  // RGT
  face = pf[3];
  reorient = get_reorient(face, cell, normalIn_H[3], M);
  children = faceTree->children(face);
  assert(children->n == 2);
  for (E_Int i = 0; i < 2; i++) RGT[i] = children->pc[i];
  Q6_get_ordered_data(M, NODES[1], reorient, RGT, local);
  assert(local[0] == NODES[1]);
  assert(local[1] == NODES[2]);
  assert(local[4] == NODES[9]);
  NODES[6] = local[2];
  NODES[5] = local[3];
  NODES[14] = local[5];

  // TOP
  face = pf[1];
  if (M->mode_2D)
    assert(check_face_aligned_with_vector(face, M->mode_2D, M));
  reorient = get_reorient(face, cell, normalIn_H[1], M);
  children = faceTree->children(face);
  assert(children->n == 4);
  for (E_Int i = 0; i < 4; i++) TOP[i] = children->pc[i];
  Q9_get_ordered_data(M, NODES[4], reorient, TOP, local);
  assert(local[0] == NODES[4]);
  assert(local[1] == NODES[5]);
  assert(local[2] == NODES[6]);
  assert(local[3] == NODES[7]);
  assert(local[5] == NODES[14]);
  assert(local[7] == NODES[13]);
  NODES[15] = local[4];
  NODES[16] = local[6];
  NODES[17] = local[8];

  // FRO
  face = pf[4];
  reorient = get_reorient(face, cell, normalIn_H[4], M);
  children = faceTree->children(face);
  assert(children->n == 2);
  for (E_Int i = 0; i < 2; i++) FRO[i] = children->pc[i];
  Q6_get_ordered_data(M, NODES[1], reorient, FRO, local);
  assert(local[0] == NODES[1]);
  assert(local[1] == NODES[0]);
  assert(local[2] == NODES[4]);
  assert(local[3] == NODES[5]);
  assert(local[4] == NODES[8]);
  assert(local[5] == NODES[15]);

  // BCK
  face = pf[5];
  reorient = get_reorient(face, cell, normalIn_H[5], M);
  children = faceTree->children(face);
  assert(children->n == 2);
  for (E_Int i = 0; i < 2; i++) BCK[i] = children->pc[i];
  Q6_get_ordered_data(M, NODES[2], reorient, BCK, local);
  assert(local[0] == NODES[2]);
  assert(local[1] == NODES[3]);
  assert(local[2] == NODES[7]);
  assert(local[3] == NODES[6]);
  assert(local[4] == NODES[10]);
  assert(local[5] == NODES[16]);

  // Set internal faces in ngon
  E_Int *ptr = &M->indPG[M->nfaces];

  // nfaces (RGT c0)
  ptr[1] = ptr[0] + 4;
  M->ngon[*ptr+0] = NODES[8];
  M->ngon[*ptr+1] = NODES[12];
  M->ngon[*ptr+2] = NODES[17];
  M->ngon[*ptr+3] = NODES[15];
  ptr++;

  // nfaces+1 (BCK c1)
  ptr[1] = ptr[0] + 4;
  M->ngon[*ptr+0] = NODES[9];
  M->ngon[*ptr+1] = NODES[12];
  M->ngon[*ptr+2] = NODES[17];
  M->ngon[*ptr+3] = NODES[14];
  ptr++;

  // nfaces+2 (LFT c2)
  ptr[1] = ptr[0] + 4;
  M->ngon[*ptr+0] = NODES[12];
  M->ngon[*ptr+1] = NODES[10];
  M->ngon[*ptr+2] = NODES[16];
  M->ngon[*ptr+3] = NODES[17];
  ptr++;

  // nfaces+3 (FRO c3)
  ptr[1] = ptr[0] + 4;
  M->ngon[*ptr+0] = NODES[12];
  M->ngon[*ptr+1] = NODES[11];
  M->ngon[*ptr+2] = NODES[13];
  M->ngon[*ptr+3] = NODES[17];

  // Assemble children

  // First child replaces cell
  ptr = &M->indPH[cell];
  M->nface[*ptr+0] = BOT[0];
  M->nface[*ptr+1] = TOP[0];
  M->nface[*ptr+2] = LFT[0];
  M->nface[*ptr+3] = M->nfaces;
  M->nface[*ptr+4] = FRO[1];
  M->nface[*ptr+5] = M->nfaces+3;

  ptr = &M->indPH[M->ncells];

  // ncells
  ptr[1] = ptr[0] + 6;
  M->nface[*ptr+0] = BOT[1];
  M->nface[*ptr+1] = TOP[1];
  M->nface[*ptr+2] = M->nfaces;
  M->nface[*ptr+3] = RGT[0];
  M->nface[*ptr+4] = FRO[0];
  M->nface[*ptr+5] = M->nfaces+1;
  ptr++;

  // ncells+1
  ptr[1] = ptr[0] + 6;
  M->nface[*ptr+0] = BOT[2];
  M->nface[*ptr+1] = TOP[2];
  M->nface[*ptr+2] = M->nfaces+2;
  M->nface[*ptr+3] = RGT[1];
  M->nface[*ptr+4] = M->nfaces+1;
  M->nface[*ptr+5] = BCK[0];
  ptr++;

  // ncells+2
  ptr[1] = ptr[0] + 6;
  M->nface[*ptr+0] = BOT[3];
  M->nface[*ptr+1] = TOP[3];
  M->nface[*ptr+2] = LFT[1];
  M->nface[*ptr+3] = M->nfaces+2;
  M->nface[*ptr+4] = M->nfaces+3;
  M->nface[*ptr+5] = BCK[1];

  // Set new cells in tree
  E_Int next_level = M->cellTree->level(cell) + 1;

  // Parent
  M->cellTree->set_parent_elem(cell, 4, M->ncells);

  // Hexa children
  for (E_Int i = 1; i < 4; i++)
    M->cellTree->set_child_elem(i, cell, HEXA, next_level, M->ncells);

  // Set external faces owns and neis
  update_external_pe_after_ref(cell, M);

  // Set internal faces in faceTree
  for (E_Int i = 0; i < 4; i++)
    M->faceTree->set_new_elem(M->nfaces+i, QUAD, next_level);

  // Set internal faces parent elements
  M->owner[M->nfaces] = cell;
  M->neigh[M->nfaces] = M->ncells;
  M->owner[M->nfaces+1] = M->ncells;
  M->neigh[M->nfaces+1] = M->ncells+1;
  M->owner[M->nfaces+2] = M->ncells+2;
  M->neigh[M->nfaces+2] = M->ncells+1;
  M->owner[M->nfaces+3] = cell;
  M->neigh[M->nfaces+3] = M->ncells+2;

  assert(check_canon_hexa(cell, M));
  for (E_Int i = 0; i < 3; i++)
    assert(check_canon_hexa(M->ncells+i, M));

  M->nfaces += 4;
  M->ncells += 3;
}

void H18_unrefine(E_Int cell, AMesh *M)
{
  assert(M->cellTree->level(cell) > 0);

  M->cellTree->state_[cell] = UNREFINED;
  M->cellTree->level_[cell] -= 1;

  Children *cur = M->cellTree->children(cell);

  assert(cur != NULL);

  for (E_Int i = 1; i < cur->n; i++) M->cellTree->state_[cur->pc[i]] = GONE;

  E_Int c0 = cur->pc[0];
  E_Int c1 = cur->pc[1];
  E_Int c2 = cur->pc[2];
  E_Int c3 = cur->pc[3];

  E_Int *f0 = get_facets(c0, M->nface, M->indPH);
  E_Int *f1 = get_facets(c1, M->nface, M->indPH);
  E_Int *f2 = get_facets(c2, M->nface, M->indPH);
  E_Int *f3 = get_facets(c3, M->nface, M->indPH);
  
  assert(f0[3] == f1[2]);
  assert(f1[5] == f2[4]);
  assert(f2[2] == f3[3]);
  assert(f3[4] == f0[5]);

  M->faceTree->state_[f0[3]] = GONE;
  M->owner[f0[3]] = -1;
  M->neigh[f0[3]] = -1;
  
  M->faceTree->state_[f1[5]] = GONE;
  M->owner[f1[5]] = -1;
  M->neigh[f1[5]] = -1;
  
  M->faceTree->state_[f2[2]] = GONE;
  M->owner[f2[2]] = -1;
  M->neigh[f2[2]] = -1;
  
  M->faceTree->state_[f3[4]] = GONE;
  M->owner[f3[4]] = -1;
  M->neigh[f3[4]] = -1;

  update_external_pe_after_unref(cell, M);

  E_Int *pf = get_facets(cell, M->nface, M->indPH);

  pf[0] = master_face(f0[0], M);
  pf[1] = master_face(f0[1], M);
  pf[2] = master_face(f0[2], M);
  pf[3] = master_face(f1[3], M);
  pf[4] = master_face(f0[4], M);
  pf[5] = master_face(f3[5], M);
}
