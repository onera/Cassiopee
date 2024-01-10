#include "Proto.h"

static
void set_new_face_in_tree(Element **tree, E_Int id, E_Int type, E_Int level)
{
  Element *elem = tree[id];
  elem->children = NULL;
  elem->nchildren = 0;
  elem->parent = -1;
  elem->position = 0;
  elem->type = type;
  elem->level = level;
}

static
void update_external_own_nei(E_Int cell, AMesh *M)
{
  Element *Elem = M->cellTree[cell];

  for (E_Int i = 0; i < Elem->nchildren; i++) {
    E_Int child = Elem->children[i];
    for (E_Int j = M->indPH[child]; j < M->indPH[child+1]; j++) {
      E_Int face = M->nface[j];
      if      (M->owner[face] == cell) M->owner[face] = child;
      else if (M->neigh[face] == cell) M->neigh[face] = child;
    }
  }
}

void get_ref_cells_and_faces(AMesh *M, std::vector<E_Int> &ref_cells,
  std::vector<E_Int> &ref_faces)
{
  ref_cells.clear();
  ref_faces.clear();

  E_Int *vcells = (E_Int *)XCALLOC(M->ncells, sizeof(E_Int));

  for (E_Int i = 0; i < M->nfaces; i++) {
    E_Int own = M->owner[i];
    if (M->ref_data[3*own] != 0) {
      ref_faces.push_back(i);
      if (!vcells[own]) {
        vcells[own] = 1;
        ref_cells.push_back(own);
      }
      continue;
    }

    E_Int nei = M->neigh[i];
    if (nei == -1) continue;

    if (M->ref_data[3*nei] != 0) {
      ref_faces.push_back(i);
      if (!vcells[nei]) {
        vcells[nei] = 1;
        ref_cells.push_back(nei);
      }
    }
  }

  XFREE(vcells);
}

static
void set_parent_elem(Element **elemTree, E_Int elem, E_Int nchildren,
  E_Int type, E_Int where)
{
  elemTree[elem]->nchildren = nchildren;

  elemTree[elem]->children = (E_Int *)XMALLOC(nchildren * sizeof(E_Int));
  elemTree[elem]->children[0] = elem;
  for (E_Int i = 1; i < nchildren; i++)
    elemTree[elem]->children[i] = where + i - 1;

  elemTree[elem]->parent = elem;
  elemTree[elem]->position = 0;
  elemTree[elem]->type = type;
  elemTree[elem]->level += 1;
}

static
void set_child_elem(Element **elemTree, E_Int parent, E_Int cpos, E_Int type,
  E_Int next_level, E_Int where)
{
  E_Int pos = where + cpos - 1;
  assert(elemTree[pos]->children == NULL);
  elemTree[pos]->nchildren = 0;
  elemTree[pos]->parent = parent;
  elemTree[pos]->position = cpos;
  elemTree[pos]->type = type;
  elemTree[pos]->level = next_level;
}

static
void refine_tri(E_Int face, AMesh *M)
{
  if (M->faceTree[face]->nchildren > 0) return;

  E_Int *pn = get_facets(face, M->ngon, M->indPG);

  Edge E;
  E_Int ec[3]; // Edge center nodes

  for (E_Int i = 0; i < 3; i++) {
    E_Int ni = pn[i];
    E_Int nj = pn[(i+1)%3];
    E.set(ni, nj);

    auto search = M->ecenter.find(E);
    if (search == M->ecenter.end()) {
      M->x[M->npoints] = 0.5*(M->x[ni] + M->x[nj]);
      M->y[M->npoints] = 0.5*(M->y[ni] + M->y[nj]);
      M->z[M->npoints] = 0.5*(M->z[ni] + M->z[nj]);
      ec[i] = M->npoints;
      M->ecenter[E] = M->npoints++;
    } else {
      ec[i] = search->second;
    }
  }

  // Set the rest of the children in ngon, starting from ptr
  E_Int *ptr = &M->indPG[M->nfaces];

  // Second face
  ptr[1] = *ptr + 3;
  M->ngon[*ptr    ] = ec[0];
  M->ngon[*ptr + 1] = pn[1];
  M->ngon[*ptr + 2] = ec[1];
  ptr++;

  // Third face
  ptr[1] = *ptr + 3;
  M->ngon[*ptr    ] = ec[2];
  M->ngon[*ptr + 1] = ec[1];
  M->ngon[*ptr + 2] = pn[2];
  ptr++;

  // Fourth face
  ptr[1] = *ptr + 3;
  M->ngon[*ptr    ] = ec[1];
  M->ngon[*ptr + 1] = ec[2];
  M->ngon[*ptr + 2] = ec[0];

  // First child replaces face
  ptr = &M->indPG[face];
  M->ngon[*ptr    ] = pn[0];
  M->ngon[*ptr + 1] = ec[0];
  M->ngon[*ptr + 2] = ec[2];

  // Set faces in faceTree
  E_Int next_level = M->faceTree[face]->level + 1;
  set_parent_elem(M->faceTree, face, 4, TRI, M->nfaces);
  for (E_Int i = 1; i < 4; i++)
    set_child_elem(M->faceTree, face, i, TRI, next_level, M->nfaces);

  // Own and nei
  M->owner[M->nfaces] = M->owner[face];
  M->neigh[M->nfaces] = M->neigh[face];
  M->owner[M->nfaces+1] = M->owner[face];
  M->neigh[M->nfaces+1] = M->neigh[face];
  M->owner[M->nfaces+2] = M->owner[face];
  M->neigh[M->nfaces+2] = M->neigh[face];

  M->nfaces += 3;
}

static
void T6_get_ordered_data(AMesh *M, E_Int NODE, E_Int reorient,
  E_Int *children, E_Int *local)
{
  E_Int *pn0 = &M->ngon[M->indPG[children[0]]];
  E_Int *pn1 = &M->ngon[M->indPG[children[1]]];
  E_Int *pn2 = &M->ngon[M->indPG[children[2]]];

  local[0] = pn0[0];
  local[1] = pn1[1];
  local[2] = pn2[2];

  local[3] = pn0[1];
  local[4] = pn1[2];
  local[5] = pn2[0];

  E_Int i0 = Get_pos(NODE, local, 3);

  assert(i0 != -1);

  Right_shift(local, i0, 3);
  Right_shift(local+3, i0, 3);
  Right_shift(children, i0, 3);

  if (reorient) {
    std::swap(children[1], children[2]);
    std::swap(local[1], local[2]);
    std::swap(local[3], local[5]);
  }
}

static
void refine_quad(E_Int face, AMesh *M)
{
  if (M->faceTree[face]->nchildren > 0) return;

  E_Int *pn = get_facets(face, M->ngon, M->indPG);

  Edge E;
  E_Int ec[4]; // Edge center nodes

  for (E_Int i = 0; i < 4; i++) {
    E_Int ni = pn[i];
    E_Int nj = pn[(i+1)%4];
    E.set(ni, nj);

    auto search = M->ecenter.find(E);
    if (search == M->ecenter.end()) {
      M->x[M->npoints] = 0.5*(M->x[ni] + M->x[nj]);
      M->y[M->npoints] = 0.5*(M->y[ni] + M->y[nj]);
      M->z[M->npoints] = 0.5*(M->z[ni] + M->z[nj]);
      ec[i] = M->npoints;
      M->ecenter[E] = M->npoints++;
    } else {
      ec[i] = search->second;
    }
  }

  // Note(Imad): face center is supposed to be already computed
  E_Float *fc = &M->fc[3*face];
  M->x[M->npoints] = fc[0];
  M->y[M->npoints] = fc[1];
  M->z[M->npoints] = fc[2];
  E_Int ncenter = M->npoints++;

  // Set the rest of the children in ngon, starting from ptr
  E_Int *ptr = &M->indPG[M->nfaces];

  // Second face
  ptr[1] = *ptr + 4;
  M->ngon[*ptr    ] = ec[0];
  M->ngon[*ptr + 1] = pn[1];
  M->ngon[*ptr + 2] = ec[1];
  M->ngon[*ptr + 3] = ncenter;
  ptr++;

  // Third face
  ptr[1] = *ptr + 4;
  M->ngon[*ptr    ] = ncenter;
  M->ngon[*ptr + 1] = ec[1];
  M->ngon[*ptr + 2] = pn[2];
  M->ngon[*ptr + 3] = ec[2];
  ptr++;

  // Fourth face
  ptr[1] = *ptr + 4;
  M->ngon[*ptr    ] = ec[3];
  M->ngon[*ptr + 1] = ncenter;
  M->ngon[*ptr + 2] = ec[2];
  M->ngon[*ptr + 3] = pn[3];

  // First child replaces face
  ptr = &M->indPG[face];
  M->ngon[*ptr    ] = pn[0];
  M->ngon[*ptr + 1] = ec[0];
  M->ngon[*ptr + 2] = ncenter;
  M->ngon[*ptr + 3] = ec[3];

  // Set faces in faceTree
  E_Int next_level = M->faceTree[face]->level + 1;
  set_parent_elem(M->faceTree, face, 4, QUAD, M->nfaces);
  for (E_Int i = 1; i < 4; i++)
    set_child_elem(M->faceTree, face, i, QUAD, next_level, M->nfaces);

  M->owner[M->nfaces] = M->owner[face];
  M->neigh[M->nfaces] = M->neigh[face];
  M->owner[M->nfaces+1] = M->owner[face];
  M->neigh[M->nfaces+1] = M->neigh[face];
  M->owner[M->nfaces+2] = M->owner[face];
  M->neigh[M->nfaces+2] = M->neigh[face];

  M->nfaces += 3;
}

static
void Q9_get_ordered_data(AMesh *M, E_Int NODE, E_Int reorient,
  E_Int *children, E_Int *local)
{
  E_Int *pn0 = &M->ngon[M->indPG[children[0]]];
  E_Int *pn1 = &M->ngon[M->indPG[children[1]]];
  E_Int *pn2 = &M->ngon[M->indPG[children[2]]];
  E_Int *pn3 = &M->ngon[M->indPG[children[3]]];

  local[0] = pn0[0];
  local[1] = pn1[1];
  local[2] = pn2[2];
  local[3] = pn3[3];

  local[4] = pn0[1];
  local[5] = pn1[2];
  local[6] = pn2[3];
  local[7] = pn3[0];

  local[8] = pn0[2];

  E_Int i0 = Get_pos(NODE, local, 4);
  assert(i0 != -1);

  Right_shift(local, i0, 4);
  Right_shift(local+4, i0, 4);
  Right_shift(children, i0, 4);

  if (reorient) {
  	std::swap(local[1], local[3]);
  	std::swap(local[4], local[7]);
  	std::swap(local[5], local[6]);
  	std::swap(children[1], children[3]);
  }
}

static inline
E_Int Tree_get_nchildren(E_Int id, Element **tree)
{
  return tree[id]->nchildren;
}

static inline
E_Int *Tree_get_children(E_Int id, Element **tree)
{
  return tree[id]->children;
}

static
void refine_hexa(E_Int cell, AMesh *M)
{
  E_Int NODES[27], FACES[24];
  for (E_Int i = 0; i < 27; i++) NODES[i] = -1;
  E_Int *BOT = FACES;
  E_Int *TOP = FACES + 4;
  E_Int *LFT = FACES + 8;
  E_Int *RGT = FACES + 12;
  E_Int *FRO = FACES + 16;
  E_Int *BCK = FACES + 20;
  E_Int *pf = get_facets(cell, M->nface, M->indPH);

  E_Int face, i0, reorient, *children, nchildren, local[9];
  Element **faceTree = M->faceTree;

  // BOT
  face = pf[0];
  reorient = get_reorient(face, cell, normalIn_H[0], M);
  children = Tree_get_children(face, M->faceTree);
  for (E_Int i = 0; i < 4; i++) BOT[i] = children[i];
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
  children = Tree_get_children(face, M->faceTree);
  for (E_Int i = 0; i < 4; i++) LFT[i] = children[i];
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
  children = Tree_get_children(face, M->faceTree);
  for (E_Int i = 0; i < 4; i++) RGT[i] = children[i];
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
  children = Tree_get_children(face, M->faceTree);
  for (E_Int i = 0; i < 4; i++) FRO[i] = children[i];
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
  children = Tree_get_children(face, M->faceTree);
  for (E_Int i = 0; i < 4; i++) BCK[i] = children[i];
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
  children = Tree_get_children(face, M->faceTree);
  for (E_Int i = 0; i < 4; i++) TOP[i] = children[i];
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

  E_Int next_level = M->cellTree[cell]->level + 1;

  // Parent
  Element *Elem = M->cellTree[cell];
  Elem->children = (E_Int *)XMALLOC(8 * sizeof(E_Int));
  Elem->children[0] = cell;
  for (E_Int i = 1; i < 8; i++) Elem->children[i] = M->ncells+i-1;
  Elem->nchildren = 8;
  Elem->parent = cell;
  Elem->position = 0;
  Elem->type = HEXA;
  Elem->level = next_level;

  // Hexa children
  for (E_Int i = 0; i < 7; i++) {
    Element *Elem = M->cellTree[M->ncells+i];
    Elem->children = NULL;
    Elem->nchildren = 0;
    Elem->parent = cell;
    Elem->position = i+1;
    Elem->type = HEXA;
    Elem->level = next_level;
  }

  //E_Int next_level = M->cellTree[cell]->level + 1;
  //set_parent_elem(M->cellTree, cell, 8, HEXA, M->ncells);
  //for (E_Int i = 0; i < 7; i++)
  //  set_child_elem(M->cellTree, cell, i, HEXA, next_level, M->ncells);

  // Set external faces owns and neis
  update_external_own_nei(cell, M);

  // Set internal faces in faceTree
  for (E_Int i = 0; i < 12; i++)
    set_new_face_in_tree(faceTree, M->nfaces+i, QUAD, next_level);

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

static
void refine_tetra(E_Int cell, AMesh *M)
{
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

  // Tetra children
  for (E_Int i = 0; i < 7; i++) {
    Element *Elem = M->cellTree[M->ncells+i];
    Elem->children = NULL;
    Elem->nchildren = 0;
    Elem->parent = cell;
    Elem->position = i+1;
    Elem->type = TETRA;
    Elem->level = next_level;
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
}

static
void refine_penta(E_Int cell, AMesh *M)
{
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

  /***************/

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

  // Penta children
  for (E_Int i = 0; i < 7; i++) {
    Element *Elem = M->cellTree[M->ncells+i];
    Elem->children = NULL;
    Elem->nchildren = 0;
    Elem->parent = cell;
    Elem->position = i+1;
    Elem->type = PENTA;
    Elem->level = next_level;
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
}

static
void refine_pyra(E_Int cell, AMesh *M)
{
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

  /*************/
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

  /*************/

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

  /*************/

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

  // Pyra children
  for (E_Int i = 0; i < 5; i++) {
    Element *Elem = M->cellTree[M->ncells+i];
    Elem->children = NULL;
    Elem->nchildren = 0;
    Elem->parent = cell;
    Elem->position = i+1;
    Elem->type = PYRA;
    Elem->level = next_level;
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
}

void refine_faces(const std::vector<E_Int> &ref_faces, AMesh *M)
{
  for (size_t i = 0; i < ref_faces.size(); i++) {
    E_Int face = ref_faces[i];
    E_Int type = M->faceTree[face]->type;
    if      (type == TRI)  refine_tri(face, M);
    else if (type == QUAD) refine_quad(face, M);
    else assert(0);
  }
}

void refine_cells(const std::vector<E_Int> &ref_cells, AMesh *M)
{
  for (size_t i = 0; i < ref_cells.size(); i++) {
    E_Int cell = ref_cells[i];
    E_Int type = M->cellTree[cell]->type;
    switch (type) {
      case HEXA:
        refine_hexa(cell, M);
        break;
      case TETRA:
        refine_tetra(cell, M);
        break;
      case PENTA:
        refine_penta(cell, M);
        break;
      case PYRA:
        refine_pyra(cell, M);
        break;
      default:
        assert(0);
        break;
    }
  }
}

void unrefine_cells(const std::vector<E_Int> &ref_cells, AMesh *M)
{
  for (size_t i = 0; i < ref_cells.size(); i++) {
    E_Int cell = ref_cells[i];
  }
}

void resize_data_for_refinement(AMesh *M, E_Int nref_cells, E_Int nref_faces)
{
  E_Int cell_incr = nref_cells * 10;
  E_Int face_incr = nref_faces * 4 // OK for quads and tris
                  + cell_incr * 13; // max 13 internal faces per refined cell
  E_Int nnew_cells = M->ncells + cell_incr;
  E_Int nnew_faces = M->nfaces + face_incr;
  // max 5 new points per refined quad + nref_cells centroids
  E_Int nnew_points = M->npoints + nref_faces*5 + nref_cells;
  
  M->cellTree = (Element **)XRESIZE(M->cellTree, nnew_cells*sizeof(Element *));
  M->faceTree = (Element **)XRESIZE(M->faceTree, nnew_faces*sizeof(Element *));
  for (E_Int i = M->ncells; i < nnew_cells; i++)
    M->cellTree[i] = (Element *)XCALLOC(1, sizeof(Element));
  for (E_Int i = M->nfaces; i < nnew_faces; i++)
    M->faceTree[i] = (Element *)XCALLOC(1, sizeof(Element));
  
  M->ngon  = (E_Int *)  XRESIZE(M->ngon,  (4*nnew_faces) * sizeof(E_Int));
  M->indPG = (E_Int *)  XRESIZE(M->indPG, (nnew_faces+1) * sizeof(E_Int));
  M->nface = (E_Int *)  XRESIZE(M->nface, (6*nnew_cells) * sizeof(E_Int));
  M->indPH = (E_Int *)  XRESIZE(M->indPH, (nnew_cells+1) * sizeof(E_Int));
  M->owner = (E_Int *)  XRESIZE(M->owner, (nnew_faces)   * sizeof(E_Int));
  M->neigh = (E_Int *)  XRESIZE(M->neigh, (nnew_faces)   * sizeof(E_Int));
  M->x     = (E_Float *)XRESIZE(M->x,     (nnew_points)  * sizeof(E_Float));
  M->y     = (E_Float *)XRESIZE(M->y,     (nnew_points)  * sizeof(E_Float));
  M->z     = (E_Float *)XRESIZE(M->z,     (nnew_points)  * sizeof(E_Float));  
  //M->fc    = (E_Float *)XRESIZE(M->fc,    (3*nnew_faces) * sizeof(E_Float));
  //M->cx    = (E_Float *)XRESIZE(M->cx,    (nnew_cells)   * sizeof(E_Float));
  //M->cy    = (E_Float *)XRESIZE(M->cy,    (nnew_cells)   * sizeof(E_Float));
  //M->cz    = (E_Float *)XRESIZE(M->cz,    (nnew_cells)   * sizeof(E_Float));

  for (E_Int i = M->nfaces; i < nnew_faces; i++) {
    M->owner[i] = -1;
    M->neigh[i] = -1;
  }
}