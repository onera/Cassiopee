#include "Proto.h"

void update_external_pe_after_unref(E_Int cell, AMesh *M)
{
  Children *children = M->cellTree->children(cell);

  for (E_Int i = 0; i < children->n; i++) {
    E_Int child = children->pc[i];

    for (E_Int j = M->indPH[child]; j < M->indPH[child+1]; j++) {
      E_Int face = M->nface[j];

      if (M->owner[face] == child) {
        M->owner[face] = cell;
      } else if (M->neigh[face] == child) {
        M->neigh[face] = cell;
      } else {
        assert(M->owner[face] == -1);
        assert(M->neigh[face] == -1);
      }
    }
  }
}

void unrefine_face(E_Int face, AMesh *M)
{
  Children *cur = M->faceTree->children(face);

  if (cur == NULL) {
    assert(M->faceTree->state(face) == GONE);
    assert(M->owner[face] == -1);
    assert(M->neigh[face] == -1);
    //printf("Skipped face %d\n", face);
    return;
  }

  assert(M->faceTree->state(face) == UNTOUCHED);

  M->faceTree->state_[face] = UNREFINED;
  M->faceTree->level_[face] -= 1;
  
  E_Int nei = M->neigh[face];
  E_Int own = M->owner[face];

  if (nei == -1) assert(master_cell(nei, M) == -1);

  if (M->cellTree->state(own) == GONE) M->owner[face] = master_cell(own, M);
  if (M->cellTree->state(nei) == GONE) M->neigh[face] = master_cell(nei, M);

  assert(M->owner[face] != M->neigh[face]);

  for (E_Int i = 1; i < cur->n; i++) {
    M->owner[cur->pc[i]] = -1;
    M->neigh[cur->pc[i]] = -1;
    M->faceTree->state_[cur->pc[i]] = GONE;
  }

  E_Int np = -1;
  E_Int *pn0 = get_face(cur->pc[0], np, M->ngon, M->indPG);
  assert(np == 4);
  E_Int *pn1 = get_face(cur->pc[1], np, M->ngon, M->indPG);
  assert(np == 4);
  E_Int *pn2 = get_face(cur->pc[2], np, M->ngon, M->indPG);
  assert(np == 4);
  E_Int *pn3 = get_face(cur->pc[3], np, M->ngon, M->indPG);
  assert(np == 4);

  pn0[1] = pn1[1];
  pn0[2] = pn2[2];
  pn0[3] = pn3[3];

  // Remove cur from tree
  /*
  Children *prev = cur->next;
  M->faceTree->children_[face] = prev;
  
  // Free
  XFREE(cur);
  */


  //M->nfaces -= 3;
}

void unrefine_hexa(E_Int cell, AMesh *M)
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

void unrefine_cells(const std::vector<E_Int> &unref_cells, size_t start,
  size_t stop, AMesh *M)
{
  for (size_t i = start; i < stop; i++) {
    E_Int cell = unref_cells[i];
    E_Int type = M->cellTree->type(cell);
    switch (type) {
      case HEXA:
        unrefine_hexa(cell, M);
        break;
      default:
        assert(0);
        break;
    }
  }
}

void unrefine_faces(const std::vector<E_Int> &unref_faces, size_t start,
  size_t stop, AMesh *M)
{
  for (size_t i = start; i < stop; i++) {
    E_Int face = unref_faces[i];
    unrefine_face(face, M);
  }
}
