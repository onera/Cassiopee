#include "proto.h"

static inline
E_Int *get_children(tree *T, E_Int id)
{
  if (T->indir[id] == -1) return NULL;
  return &T->children[T->indir[id]+1];
}

static inline
E_Int get_nchildren(tree *T, E_Int id)
{
  if (T->indir[id] == -1) return 0;
  return T->children[T->indir[id]];
}

void resize_data_for_refinement(mesh *M, tree *ct, tree *ft, E_Int nref_cells)
{
  E_Int cell_increment = 8*nref_cells; // 8 new cells per cell
  E_Int face_increment = 36*nref_cells; // 24 cfaces + 12 ifaces per new cell

  tree_resize(ft, face_increment, 4);
  tree_resize(ct, cell_increment, 8);

  E_Int new_nfaces = M->nfaces + face_increment;
  E_Int new_ncells = M->ncells + cell_increment;
	E_Int new_npoints = M->npoints + cell_increment*19; // 27 - 8

  M->NGON  = (E_Int *)   XRESIZE(M->NGON,   (4*new_nfaces)  * sizeof(E_Int));
  M->NFACE  = (E_Int *)  XRESIZE(M->NFACE,  (6*new_ncells)  * sizeof(E_Int));
  M->fc     = (E_Float *)XRESIZE(M->fc,     (3*new_nfaces)  * sizeof(E_Float));
  M->cc     = (E_Float *)XRESIZE(M->cc,     (3*new_ncells)  * sizeof(E_Float));
  M->xyz    = (E_Float *)XRESIZE(M->xyz,    (3*new_npoints) * sizeof(E_Float));
  M->owner  = (E_Int *)  XRESIZE(M->owner,  new_nfaces      * sizeof(E_Int));
  M->neigh  = (E_Int *)  XRESIZE(M->neigh,  new_nfaces      * sizeof(E_Int));
  M->xfaces = (E_Int *)  XRESIZE(M->xfaces, (new_nfaces+1)  * sizeof(E_Int));
  M->xcells = (E_Int *)  XRESIZE(M->xcells, (new_ncells+1)  * sizeof(E_Int));

  M->predicted_ncells = M->ncells;
  M->predicted_nfaces = M->nfaces;
  M->predicted_npoints = M->npoints;

  // init
  for (E_Int i = M->nfaces; i < new_nfaces; i++) {
    M->owner[i] = -1;
    M->neigh[i] = -1;
  }

  // assume quad and hexa
  for (E_Int i = M->nfaces; i < new_nfaces; i++)
    M->xfaces[i+1] = 4;
  for (E_Int i = M->nfaces; i < new_nfaces; i++)
    M->xfaces[i+1] += M->xfaces[i];

  for (E_Int i = M->ncells; i < new_ncells; i++)
    M->xcells[i+1] = 6;
  for (E_Int i = M->ncells; i < new_ncells; i++)
    M->xcells[i+1] += M->xcells[i];
}

E_Int is_cell_to_refine(E_Int *ref_data)
{
  return (ref_data[0] > 0) || (ref_data[1] > 0) || (ref_data[2] > 0);
}

std::vector<E_Int> get_ref_cells(mesh *M)
{
  std::vector<E_Int> ref_cells;
  for (E_Int i = 0; i < M->ncells; i++) {
    if (is_cell_to_refine(&M->ref_data[3*i]))
      ref_cells.push_back(i);
  }

  return ref_cells;
}

void cut_face_x(E_Int face, mesh *M, tree *ft, E_Int i0, E_Int reorient)
{
  assert(get_nchildren(ft, face) == 0 || get_nchildren(ft, face) == 2);
  if (get_nchildren(ft, face) == 2) return;
  
  if (reorient) {
    if (i0 == 0 || i0 == 2) {
      cut_face_y(face, M, ft, 0, 0);
      return;
    }
  } else {
    if (i0 == 1 || i0 == 3) {
      cut_face_y(face, M, ft, 0, 0);
      return;
    }
  }

  E_Int i;
  E_Int *faces = M->NGON;
  E_Int *pn = &M->NGON[4*face];
  E_Int n0 = pn[0];
  E_Int n1 = pn[1];
  E_Int n2 = pn[2];
  E_Int n3 = pn[3];

  E_Float *xyz = M->xyz;
  E_Float *pt, *x0, *x1;
  E_Int n4, n5;
  
  edge E;

  E.set(n0, n1);
  {
    auto search = M->ET.find(E);
    if (search == M->ET.end()) {
      pt = &xyz[3*M->npoints];
      x0 = &xyz[3*n0];
      x1 = &xyz[3*n1];
      for (i = 0; i < 3; i++)
        pt[i] = 0.5*(x0[i] + x1[i]);
      n4 = M->npoints;
      M->ET[E] = M->npoints++;
    } else {
      n4 = search->second;
    }
  }

  E.set(n2, n3);
  {
    auto search = M->ET.find(E);
    if (search == M->ET.end()) {
      pt = &xyz[3*M->npoints];
      x0 = &xyz[3*n2];
      x1 = &xyz[3*n3];
      for (i = 0; i < 3; i++)
        pt[i] = 0.5*(x0[i] + x1[i]);
      n5 = M->npoints;
      M->ET[E] = M->npoints++;
    } else {
      n5 = search->second;
    }
  }

  E_Int *new_face_0 = &faces[4*M->nfaces];
  E_Int *new_face_1 = &faces[4*(M->nfaces+1)];
  
  new_face_0[0] = n0; new_face_1[0] = n4;
  new_face_0[1] = n4; new_face_1[1] = n1;
  new_face_0[2] = n5; new_face_1[2] = n2;
  new_face_0[3] = n3; new_face_1[3] = n5;
  
  // compute new face centers
  for (i = 0; i < 2; i++)
  	compute_face_center(M, M->nfaces+i);
  
  tree_insert_children(ft, face, M->nfaces, 2);
  
  M->owner[M->nfaces] = M->owner[face];
  M->neigh[M->nfaces] = M->neigh[face];
  M->owner[M->nfaces+1] = M->owner[face];
  M->neigh[M->nfaces+1] = M->neigh[face];
  M->nfaces += 2;
}

void cut_face_y(E_Int face, mesh *M, tree *ft, E_Int i0, E_Int reorient)
{
  assert(get_nchildren(ft, face) == 0 || get_nchildren(ft, face) == 2);
  
  if (get_nchildren(ft, face) == 2) return; 
  
  if (reorient) {
    if (i0 == 0 || i0 == 2) {
  	  cut_face_x(face, M, ft, 0, 0);
  	  return;
    }
  } else {
    if (i0 == 1 || i0 == 3) {
      cut_face_x(face, M, ft, 0, 0);
      return;
    }
  } 
  
  E_Int i, *faces;
  faces = M->NGON;
  E_Int *pn = &faces[4*face];
  E_Float *xyz = M->xyz;  
  E_Int n0 = pn[0];
  E_Int n1 = pn[1];
  E_Int n2 = pn[2];
  E_Int n3 = pn[3]; 
  E_Float *pt, *x0, *x1;
  E_Int n4, n5; 
  edge E;

  E.set(n1, n2);
  {
    auto search = M->ET.find(E);
    if (search == M->ET.end()) {
      pt = &xyz[3*M->npoints];
      x0 = &xyz[3*n1];
      x1 = &xyz[3*n2];
      for (i = 0; i < 3; i++)
        pt[i] = 0.5*(x0[i] + x1[i]);
      n4 = M->npoints;
      M->ET[E] = M->npoints++;
    } else {
      n4 = search->second;
    }
  }

  E.set(n3, n0);
  {
    auto search = M->ET.find(E);
    if (search == M->ET.end()) {
      pt = &xyz[3*M->npoints];
      x0 = &xyz[3*n3];
      x1 = &xyz[3*n0];
      for (i = 0; i < 3; i++)
        pt[i] = 0.5*(x0[i] + x1[i]);
      n5 = M->npoints;
      M->ET[E] = M->npoints++;
    } else {
      n5 = search->second;
    }
  }

  E_Int *new_face_0 = &M->NGON[4*M->nfaces];
  E_Int *new_face_1 = &M->NGON[4*(M->nfaces+1)];
  
  new_face_0[0] = n0; new_face_1[0] = n5;
  new_face_0[1] = n1; new_face_1[1] = n4;
  new_face_0[2] = n4; new_face_1[2] = n2;
  new_face_0[3] = n5; new_face_1[3] = n3;
  
  // compute new face centers
  for (i = 0; i < 2; i++)
  	compute_face_center(M, M->nfaces+i);
  
  tree_insert_children(ft, face, M->nfaces, 2);
  
  M->owner[M->nfaces] = M->owner[face];
  M->neigh[M->nfaces] = M->neigh[face];
  M->owner[M->nfaces+1] = M->owner[face];
  M->neigh[M->nfaces+1] = M->neigh[face];
  
  M->nfaces += 2;
}

void cut_face_xy(E_Int face, mesh *M, tree *ft)
{
  assert(get_nchildren(ft, face) == 0 || get_nchildren(ft, face) == 4);

  if (get_nchildren(ft, face) == 4) return;

  E_Int i, *faces;
  faces = M->NGON;
  E_Int *pn = &faces[4*face];
  E_Float *xyz = M->xyz;
  E_Float *pt, *x0, *x1;
  E_Int n0, n1, n2, n3, n4, n[4];
  edge E;
  
  for (i = 0; i < 4; i++) {
    n0 = pn[i];
    n1 = pn[(i+1)%4];
    E.set(n0, n1);
    auto search = M->ET.find(E);
    if (search == M->ET.end()) {
      pt = &xyz[3*M->npoints];
      x0 = &xyz[3*n0];
      x1 = &xyz[3*n1];
      for (E_Int k = 0; k < 3; k++)
      	pt[k] = 0.5*(x0[k] + x1[k]);
      n[i] = M->npoints;
      M->ET[E] = M->npoints++;
    } else {
      n[i] = search->second;
    }
  }

  // NOTE: face centroid is supposed to be already computed
  E_Float *fx = &M->fc[3*face];
  E_Float *X = &xyz[3*M->npoints];
  for (i = 0; i < 3; i++)
  	X[i] = fx[i];
  
  n4 = M->npoints++;
  
  E_Int *new_face_0 = &faces[4*M->nfaces];
  E_Int *new_face_1 = &faces[4*(M->nfaces+1)];
  E_Int *new_face_2 = &faces[4*(M->nfaces+2)];
  E_Int *new_face_3 = &faces[4*(M->nfaces+3)];
  
  n0 = pn[0]; n1 = pn[1]; n2 = pn[2]; n3 = pn[3];
  
  new_face_0[0] = n0;   new_face_1[0] = n[0];
  new_face_0[1] = n[0]; new_face_1[1] = n1;
  new_face_0[2] = n4;   new_face_1[2] = n[1];
  new_face_0[3] = n[3]; new_face_1[3] = n4;
  new_face_2[0] = n4;   new_face_3[0] = n[3];
  new_face_2[1] = n[1]; new_face_3[1] = n4;
  new_face_2[2] = n2;   new_face_3[2] = n[2];
  new_face_2[3] = n[2]; new_face_3[3] = n3;
  
  // compute new face centers
  for (i = 0; i < 4; i++)
  	compute_face_center(M, M->nfaces+i);
  
  tree_insert_children(ft, face, M->nfaces, 4);
  
  E_Int own = M->owner[face];
  E_Int nei = M->neigh[face];
  
  M->owner[M->nfaces] = own;
  M->neigh[M->nfaces] = nei;
  M->owner[M->nfaces+1] = own;
  M->neigh[M->nfaces+1] = nei;
  M->owner[M->nfaces+2] = own;
  M->neigh[M->nfaces+2] = nei;
  M->owner[M->nfaces+3] = own;
  M->neigh[M->nfaces+3] = nei;
  
  M->nfaces += 4;
}

static
void update_external_own_nei(E_Int cell, mesh *M, tree *ct, tree *ft)
{
  E_Int i, j, child, *pf, face, *cells;
  cells = M->NFACE;
  E_Int *children = get_children(ct, cell);
  E_Int nchildren = get_nchildren(ct, cell);
  for (i = 0; i < nchildren; i++) {
    child = children[i];
    pf = &cells[6*child];
    for (j = 0; j < 6; j++) {
      face = pf[j];
      if (M->owner[face] == -1 && M->neigh[face] == -1)
        continue;
      if (M->owner[face] == cell) M->owner[face] = child;
      else M->neigh[face] = child;
    }
  }
}

static
void Q6_get_ordered_data(E_Int face, mesh *M, E_Int i0, E_Int reorient,
  E_Int *children, E_Int *local)
{
  E_Int *pn = &M->NGON[4*face];
  for (E_Int i = 0; i < 4; i++) local[i] = pn[i];
  E_Int *pn0 = &M->NGON[4*children[0]];
  E_Int dir = (pn0[1] == pn[1]) ? 1 : 0; // X = 0, Y = 1
  right_shift(local, i0, 4);
  if (reorient) std::swap(local[1], local[3]);

  if (dir == 0) {
    if (!reorient) {
      switch (i0) {
        case 0:
        	local[4] = pn0[1];
          local[5] = pn0[2];
        	break;
        case 1:
          local[4] = pn0[2];
          local[5] = pn0[1];
          std::swap(children[0], children[1]);
        	break;
        case 2:
          local[4] = pn0[2];
          local[5] = pn0[1];
          std::swap(children[0], children[1]);
        	break;
        case 3:
        	local[4] = pn0[1];
          local[5] = pn0[2];
        	break;
      }
    } else {
      switch (i0) {
        case 0:
        	local[4] = pn0[2];
          local[5] = pn0[1];
        	break;
        case 1:
          local[4] = pn0[1];
          local[5] = pn0[2];
          std::swap(children[0], children[1]);
        	break;
        case 2:
          local[4] = pn0[1];
          local[5] = pn0[2];
          std::swap(children[0], children[1]);
        	break;
        case 3:
        	local[4] = pn0[2];
          local[5] = pn0[1];
        	break;
      }
    }
  } else { // dir = 1
    if (!reorient) {
      switch (i0) {
        case 0:
          local[4] = pn0[2];
          local[5] = pn0[3];
          break;
        case 1:
          local[4] = pn0[2];
          local[5] = pn0[3];
          break;
        case 2:
          local[4] = pn0[3];
          local[5] = pn0[2];
          std::swap(children[0], children[1]);
          break;
        case 3:
          local[4] = pn0[3];
          local[5] = pn0[2];
          std::swap(children[0], children[1]);
          break;
      }
    } else {
      switch (i0) {
        case 0:
          local[4] = pn0[3];
          local[5] = pn0[2];
          break;
        case 1:
          local[4] = pn0[3];
          local[5] = pn0[2];
          break;
        case 2:
          local[4] = pn0[2];
          local[5] = pn0[3];
          std::swap(children[0], children[1]);
          break;
        case 3:
          local[4] = pn0[2];
          local[5] = pn0[3];
          std::swap(children[0], children[1]);
          break;
      }
    }
  }
}

static
void Q9_get_ordered_data(E_Int face, mesh *M, E_Int i0, E_Int reorient, E_Int *children, E_Int *local)
{
  E_Int *pn = &M->NGON[4*face];
  E_Int *pn0 = &M->NGON[4*children[0]];
  E_Int *pn2 = &M->NGON[4*children[2]];

  for (E_Int i = 0; i < 4; i++) local[i] = pn[i];

  local[4] = pn0[1];
  local[5] = pn2[1];
  local[6] = pn2[3];
  local[7] = pn0[3];
  local[8] = pn0[2];

  right_shift(&local[0], i0, 4);
  right_shift(&local[4], i0, 4);
  right_shift(&children[0], i0, 4);

  if (reorient) {
  	std::swap(local[1], local[3]);
  	std::swap(local[4], local[7]);
  	std::swap(local[5], local[6]);
  	std::swap(children[1], children[3]);
  }
}

void cut_cell_x(E_Int cell, mesh *M, tree *ct, tree *ft)
{
  E_Int i, face, *pn, i0, reorient, *children, local[6];
  E_Int NODES[12] = {-1};
  E_Int FACES[10];
  E_Int *BOT = FACES;
  E_Int *TOP = FACES + 2;
  E_Int *LFT = FACES + 4;
  E_Int *RGT = FACES + 5;
  E_Int *FRO = FACES + 6;
  E_Int *BCK = FACES + 8;
  E_Int *pf = &M->NFACE[6*cell];
  E_Int *faces = M->NGON;
  
  // BOT
  face = pf[0];
  pn = &faces[4*face];
  i0 = 0;
  reorient = get_reorient(face, cell, normalIn[0], M);
  cut_face_x(face, M, ft, i0, reorient);
  assert(get_nchildren(ft, face) == 2);
  children = get_children(ft, face);
  for (i = 0; i < 2; i++) BOT[i] = children[i];
  Q6_get_ordered_data(face, M, i0, reorient, BOT, local);
  for (i = 0; i < 4; i++) NODES[i] = local[i];
  NODES[8] = local[4];
  NODES[9] = local[5];
  
  // FRO
  face = pf[4];
  pn = &faces[4*face];
  i0 = get_pos(NODES[1], pn, 4);
  reorient = get_reorient(face, cell, normalIn[4], M);
  cut_face_x(face, M, ft, i0, reorient);
  assert(get_nchildren(ft, face) == 2);
  children = get_children(ft, face);
  for (i = 0; i < 2; i++) FRO[i] = children[i];
  Q6_get_ordered_data(face, M, i0, reorient, FRO, local);
  assert(local[0] == NODES[1]);
  assert(local[1] == NODES[0]);
  assert(local[4] == NODES[8]);
  NODES[4] = local[2];
  NODES[5] = local[3];
  NODES[10] = local[5];
  
  // BCK
  face = pf[5];
  pn = &faces[4*face];
  i0 = get_pos(NODES[2], pn, 4);
  reorient = get_reorient(face, cell, normalIn[5], M);
  cut_face_x(face, M, ft, i0, reorient);
  assert(get_nchildren(ft, face) == 2);
  children = get_children(ft, face);
  for (i = 0; i < 2; i++) BCK[i] = children[i];
  Q6_get_ordered_data(face, M, i0, reorient, BCK, local);
  assert(local[0] == NODES[2]);
  assert(local[1] == NODES[3]);
  assert(local[4] == NODES[9]);
  NODES[7] = local[2];
  NODES[6] = local[3];
  NODES[11] = local[5];
  
  // TOP
  face = pf[1];
  pn = &faces[4*face];
  i0 = get_pos(NODES[4], pn, 4);
  reorient = get_reorient(face, cell, normalIn[1], M);
  cut_face_x(face, M, ft, i0, reorient);
  assert(get_nchildren(ft, face) == 2);
  children = get_children(ft, face);
  for (i = 0; i < 2; i++) TOP[i] = children[i];
  Q6_get_ordered_data(face, M, i0, reorient, TOP, local);
  assert(local[0] == NODES[4]);
  assert(local[1] == NODES[5]);
  assert(local[2] == NODES[6]);
  assert(local[3] == NODES[7]);
  assert(local[4] == NODES[10]);
  assert(local[5] == NODES[11]);
  
  // LFT
  *LFT = pf[2];
  
  // RGT
  *RGT = pf[3];
  
  // E_Internal face
  E_Int nfaces = M->nfaces;
  E_Int *new_face = &faces[4*nfaces];
  new_face[0] = NODES[8];
  new_face[1] = NODES[9];
  new_face[2] = NODES[11];
  new_face[3] = NODES[10];
  compute_face_center(M, M->nfaces);
  
  // assemble children
  E_Int ncells = M->ncells;
  E_Int *new_cell_0 = &M->NFACE[6*ncells];
  E_Int *new_cell_1 = &M->NFACE[6*(ncells+1)];
  new_cell_0[0] = BOT[0]; new_cell_1[0] = BOT[1];
  new_cell_0[1] = TOP[0]; new_cell_1[1] = TOP[1];
  new_cell_0[2] = LFT[0]; new_cell_1[2] = nfaces;
  new_cell_0[3] = nfaces; new_cell_1[3] = RGT[0];
  new_cell_0[4] = FRO[1]; new_cell_1[4] = FRO[0];
  new_cell_0[5] = BCK[1]; new_cell_1[5] = BCK[0];
  
  // compute new cell centers
  for (i = 0; i < 2; i++)
  	compute_cell_center(M, M->ncells+i);
  
  // update cell tree
  tree_insert_children(ct, cell, M->ncells, 2);
  
  update_external_own_nei(cell, M, ct, ft);
  
  // set E_Internal face M->owner and M->neighbour
  M->owner[nfaces] = ncells;
  M->neigh[nfaces] = ncells+1;
  M->ncells += 2;
  M->nfaces += 1;
}

void cut_cell_y(E_Int cell, mesh *M, tree *ct, tree *ft)
{
  E_Int i, face, *pn, i0, reorient, *children, local[6];
  E_Int NODES[12] = {-1};
  E_Int FACES[10];
  E_Int *BOT = FACES;
  E_Int *TOP = FACES + 2;
  E_Int *LFT = FACES + 4;
  E_Int *RGT = FACES + 6;
  E_Int *FRO = FACES + 8;
  E_Int *BCK = FACES + 9;
  E_Int *pf = &M->NFACE[6*cell];
  E_Int *faces = M->NGON;
  
  // BOT
  face = pf[0];
  pn = &faces[4*face];
  i0 = 0;
  reorient = get_reorient(face, cell, normalIn[0], M);
  cut_face_y(face, M, ft, i0, reorient);
  assert(get_nchildren(ft, face) == 2);
  children = get_children(ft, face);
  for (i = 0; i < 2; i++) BOT[i] = children[i];
  Q6_get_ordered_data(face, M, i0, reorient, BOT, local);
  for (i = 0; i < 4; i++) NODES[i] = local[i];
  NODES[8] = local[4];
  NODES[9] = local[5];
  
  // LFT
  face = pf[2];
  pn = &faces[4*face];
  i0 = get_pos(NODES[0], pn, 4);
  reorient = get_reorient(face, cell, normalIn[2], M);
  cut_face_x(face, M, ft, i0, reorient);
  assert(get_nchildren(ft, face) == 2);
  children = get_children(ft, face);
  for (i = 0; i < 2; i++) LFT[i] = children[i];
  Q6_get_ordered_data(face, M, i0, reorient, LFT, local);
  assert(local[0] == NODES[0]);
  assert(local[1] == NODES[3]);
  assert(local[4] == NODES[9]);
  NODES[7] = local[2];
  NODES[4] = local[3];
  NODES[11] = local[5];
  
  // RGT
  face = pf[3];
  pn = &faces[4*face];
  i0 = get_pos(NODES[1], pn, 4);
  reorient = get_reorient(face, cell, normalIn[3], M);
  cut_face_x(face, M, ft, i0, reorient);
  assert(get_nchildren(ft, face) == 2);
  children = get_children(ft, face);
  for (i = 0; i < 2; i++) RGT[i] = children[i];
  Q6_get_ordered_data(face, M, i0, reorient, RGT, local);
  assert(local[0] == NODES[1]);
  assert(local[1] == NODES[2]);
  assert(local[4] == NODES[8]);
  NODES[6] = local[2];
  NODES[5] = local[3];
  NODES[10] = local[5];
  
  // TOP
  face = pf[1];
  pn = &faces[4*face];
  i0 = get_pos(NODES[4], pn, 4);
  reorient = get_reorient(face, cell, normalIn[1], M);
  cut_face_y(face, M, ft, i0, reorient);
  assert(get_nchildren(ft, face) == 2);
  children = get_children(ft, face);
  for (i = 0; i < 2; i++) TOP[i] = children[i];
  Q6_get_ordered_data(face, M, i0, reorient, TOP, local);
  assert(local[0] == NODES[4]);
  assert(local[1] == NODES[5]);
  assert(local[2] == NODES[6]);
  assert(local[3] == NODES[7]);
  assert(local[4] == NODES[10]);
  assert(local[5] == NODES[11]);
  
  // FRO
  *FRO = pf[4];
  
  // BCK
  *BCK = pf[5];
  
  // E_Internal face
  E_Int nfaces = M->nfaces;
  E_Int *new_face = &faces[4*nfaces];
  new_face[0] = NODES[8];
  new_face[1] = NODES[9];
  new_face[2] = NODES[11];
  new_face[3] = NODES[10];
  compute_face_center(M, M->nfaces);
  
  // assemble children
  E_Int ncells = M->ncells;
  E_Int *new_cell_0 = &M->NFACE[6*ncells];
  E_Int *new_cell_1 = &M->NFACE[6*(ncells+1)];
  new_cell_0[0] = BOT[0]; new_cell_1[0] = BOT[1];
  new_cell_0[1] = TOP[0]; new_cell_1[1] = TOP[1];
  new_cell_0[2] = LFT[0]; new_cell_1[2] = LFT[1];
  new_cell_0[3] = RGT[0]; new_cell_1[3] = RGT[1];
  new_cell_0[4] = FRO[0]; new_cell_1[4] = nfaces;
  new_cell_0[5] = nfaces; new_cell_1[5] = BCK[0];
  
  // compute new cell centers
  for (i = 0; i < 2; i++)
  	compute_cell_center(M, M->ncells+i);
  
  // update cell tree
  tree_insert_children(ct, cell, M->ncells, 2);
  update_external_own_nei(cell, M, ct, ft);
  
  // set E_Internal face M->owner and M->neighbour
  M->owner[nfaces] = ncells;
  M->neigh[nfaces] = ncells+1;
  M->ncells += 2;
  M->nfaces += 1;
}

void cut_cell_z(E_Int cell, mesh *M, tree *ct, tree *ft)
{
  E_Int i, face, *pn, i0, reorient, *children, local[6];
  E_Int NODES[12] = {-1};
  E_Int FACES[10];
  E_Int *BOT = FACES;
  E_Int *TOP = FACES + 1;
  E_Int *LFT = FACES + 2;
  E_Int *RGT = FACES + 4;
  E_Int *FRO = FACES + 6;
  E_Int *BCK = FACES + 8;
  E_Int *pf = &M->NFACE[6*cell];
  E_Int *faces = M->NGON;

  // BOT
  face = pf[0];
  *BOT = face;
  pn = &faces[4*face];
  i0 = 0;
  reorient = get_reorient(face, cell, normalIn[0], M);
  order_quad(local, pn, reorient, i0);
  for (i = 0; i < 4; i++) NODES[i] = local[i];
  
  // LFT
  face = pf[2];
  pn = &faces[4*face];
  i0 = get_pos(NODES[0], pn, 4);
  reorient = get_reorient(face, cell, normalIn[2], M);
  cut_face_y(face, M, ft, i0, reorient);
  assert(get_nchildren(ft, face) == 2);
  children = get_children(ft, face);
  for (i = 0; i < 2; i++) LFT[i] = children[i];
  Q6_get_ordered_data(face, M, i0, reorient, LFT, local);
  assert(local[0] == NODES[0]);
  assert(local[1] == NODES[3]);
  NODES[7] = local[2];
  NODES[4] = local[3];
  NODES[11] = local[4];
  NODES[8] = local[5];
  
  // RGT
  face = pf[3];
  pn = &faces[4*face];
  i0 = get_pos(NODES[1], pn, 4);
  reorient = get_reorient(face, cell, normalIn[3], M);
  cut_face_y(face, M, ft, i0, reorient);
  assert(get_nchildren(ft, face) == 2);
  children = get_children(ft, face);
  for (i = 0; i < 2; i++) RGT[i] = children[i];
  Q6_get_ordered_data(face, M, i0, reorient, RGT, local);
  assert(local[0] == NODES[1]);
  assert(local[1] == NODES[2]);
  NODES[6] = local[2];
  NODES[5] = local[3];
  NODES[10] = local[4];
  NODES[9] = local[5];
  
  // FRO
  face = pf[4];
  pn = &faces[4*face];
  i0 = get_pos(NODES[1], pn, 4);
  reorient = get_reorient(face, cell, normalIn[4], M);
  cut_face_y(face, M, ft, i0, reorient);
  assert(get_nchildren(ft, face) == 2);
  children = get_children(ft, face);
  for (i = 0; i < 2; i++) FRO[i] = children[i];
  Q6_get_ordered_data(face, M, i0, reorient, FRO, local);
  assert(local[0] == NODES[1]);
  assert(local[1] == NODES[0]);
  assert(local[2] == NODES[4]);
  assert(local[3] == NODES[5]);
  assert(local[4] == NODES[8]);
  assert(local[5] == NODES[9]);
  
  // BCK
  face = pf[5];
  pn = &faces[4*face];
  i0 = get_pos(NODES[2], pn, 4);
  reorient = get_reorient(face, cell, normalIn[5], M);
  cut_face_y(face, M, ft, i0, reorient);
  assert(get_nchildren(ft, face) == 2);
  children = get_children(ft, face);
  for (i = 0; i < 2; i++) BCK[i] = children[i];
  Q6_get_ordered_data(face, M, i0, reorient, BCK, local);
  assert(local[0] == NODES[2]);
  assert(local[1] == NODES[3]);
  assert(local[2] == NODES[7]);
  assert(local[3] == NODES[6]);
  assert(local[4] == NODES[11]);
  assert(local[5] == NODES[10]);
  
  // TOP
  *TOP = pf[1];
  
  // E_Internal face
  E_Int nfaces = M->nfaces;
  E_Int *new_face = &faces[4*nfaces];
  new_face[0] = NODES[8];
  new_face[1] = NODES[9];
  new_face[2] = NODES[10];
  new_face[3] = NODES[11];
  compute_face_center(M, M->nfaces);
  
  // assemble children
  E_Int ncells = M->ncells;
  E_Int *new_cell_0 = &M->NFACE[6*ncells];
  E_Int *new_cell_1 = &M->NFACE[6*(ncells+1)];
  new_cell_0[0] = BOT[0]; new_cell_1[0] = nfaces;
  new_cell_0[1] = nfaces; new_cell_1[1] = TOP[0];
  new_cell_0[2] = LFT[0]; new_cell_1[2] = LFT[1];
  new_cell_0[3] = RGT[0]; new_cell_1[3] = RGT[1];
  new_cell_0[4] = FRO[0]; new_cell_1[4] = FRO[1];
  new_cell_0[5] = BCK[0]; new_cell_1[5] = BCK[1];
  
  // compute new cell centers
  for (i = 0; i < 2; i++)
  	compute_cell_center(M, M->ncells+i);
  
  // update cell tree
  tree_insert_children(ct, cell, M->ncells, 2);
  
  update_external_own_nei(cell, M, ct, ft);
  
  // set E_Internal face M->owner and M->neighbour
  M->owner[nfaces] = ncells;
  M->neigh[nfaces] = ncells+1;
  M->ncells += 2;
  M->nfaces += 1;
}

void cut_cell_xy(E_Int cell, mesh *M, tree *ct, tree *ft)
{
  E_Int i, face, *pn, i0, reorient, *children, local[9];
  E_Int NODES[18] = {-1};
  E_Int FACES[16];
  E_Int *BOT = FACES;
  E_Int *TOP = FACES + 4;
  E_Int *LFT = FACES + 8;
  E_Int *RGT = FACES + 10;
  E_Int *FRO = FACES + 12;
  E_Int *BCK = FACES + 14;
  E_Int *pf = &M->NFACE[6*cell];
  E_Int *faces = M->NGON;
  
  // BOT
  face = pf[0];
  pn = &faces[4*face];
  i0 = 0;
  reorient = get_reorient(face, cell, normalIn[0], M);
  cut_face_xy(face, M, ft);
  assert(get_nchildren(ft, face) == 4);
  children = get_children(ft, face);
  for (i = 0; i < 4; i++) BOT[i] = children[i];
  Q9_get_ordered_data(face, M, i0, reorient, BOT, local);
  for (i = 0; i < 4; i++) NODES[i] = local[i];
  NODES[8] = local[4];
  NODES[9] = local[5];
  NODES[10] = local[6];
  NODES[11] = local[7];
  NODES[12] = local[8];
  
  // LFT
  face = pf[2];
  pn = &faces[4*face];
  i0 = get_pos(NODES[0], pn, 4);
  reorient = get_reorient(face, cell, normalIn[2], M);
  cut_face_x(face, M, ft, i0, reorient);
  assert(get_nchildren(ft, face) == 2);
  children = get_children(ft, face);
  for (i = 0; i < 2; i++) LFT[i] = children[i];
  Q6_get_ordered_data(face, M, i0, reorient, LFT, local);
  assert(local[0] == NODES[0]);
  assert(local[1] == NODES[3]);
  assert(local[4] == NODES[11]);
  NODES[7] = local[2];
  NODES[4] = local[3];
  NODES[16] = local[5];
  
  // RGT
  face = pf[3];
  pn = &faces[4*face];
  i0 = get_pos(NODES[1], pn, 4);
  reorient = get_reorient(face, cell, normalIn[3], M);
  cut_face_x(face, M, ft, i0, reorient);
  assert(get_nchildren(ft, face) == 2);
  children = get_children(ft, face);
  for (i = 0; i < 2; i++) RGT[i] = children[i];
  Q6_get_ordered_data(face, M, i0, reorient, RGT, local);
  assert(local[0] == NODES[1]);
  assert(local[1] == NODES[2]);
  assert(local[4] == NODES[9]);
  NODES[6] = local[2];
  NODES[5] = local[3];
  NODES[14] = local[5];
  
  // TOP
  face = pf[1];
  pn = &faces[4*face];
  i0 = get_pos(NODES[4], pn, 4);
  reorient = get_reorient(face, cell, normalIn[1], M);
  cut_face_xy(face, M, ft);
  assert(get_nchildren(ft, face) == 4);
  children = get_children(ft, face);
  for (i = 0; i < 4; i++) TOP[i] = children[i];
  Q9_get_ordered_data(face, M, i0, reorient, TOP, local);
  assert(local[0] == NODES[4]);
  assert(local[1] == NODES[5]);
  assert(local[2] == NODES[6]);
  assert(local[3] == NODES[7]);
  assert(local[5] == NODES[14]);
  assert(local[7] == NODES[16]);
  NODES[13] = local[4];
  NODES[15] = local[6];
  NODES[17] = local[8];
  
  // FRO
  face = pf[4];
  pn = &faces[4*face];
  i0 = get_pos(NODES[1], pn, 4);
  reorient = get_reorient(face, cell, normalIn[4], M);
  cut_face_x(face, M, ft, i0, reorient);
  assert(get_nchildren(ft, face) == 2);
  children = get_children(ft, face);
  for (i = 0; i < 2; i++) FRO[i] = children[i];
  Q6_get_ordered_data(face, M, i0, reorient, FRO, local);
  assert(local[0] == NODES[1]);
  assert(local[1] == NODES[0]);
  assert(local[2] == NODES[4]);
  assert(local[3] == NODES[5]);
  assert(local[4] == NODES[8]);
  assert(local[5] == NODES[13]);
  
  // BCK
  face = pf[5];
  pn = &faces[4*face];
  i0 = get_pos(NODES[2], pn, 4);
  reorient = get_reorient(face, cell, normalIn[5], M);
  cut_face_x(face, M, ft, i0, reorient);
  assert(get_nchildren(ft, face) == 2);
  children = get_children(ft, face);
  for (i = 0; i < 2; i++) BCK[i] = children[i];
  Q6_get_ordered_data(face, M, i0, reorient, BCK, local);
  assert(local[0] == NODES[2]);
  assert(local[1] == NODES[3]);
  assert(local[2] == NODES[7]);
  assert(local[3] == NODES[6]);
  assert(local[4] == NODES[10]);
  assert(local[5] == NODES[15]);
  
  // E_Internal faces
  E_Int nfaces = M->nfaces;
  E_Int *new_face_0 = &faces[4*nfaces];
  E_Int *new_face_1 = &faces[4*(nfaces+1)];
  E_Int *new_face_2 = &faces[4*(nfaces+2)];
  E_Int *new_face_3 = &faces[4*(nfaces+3)];
  new_face_0[0] = NODES[8];   new_face_1[0] = NODES[12];
  new_face_0[1] = NODES[12];  new_face_1[1] = NODES[10];
  new_face_0[2] = NODES[17];  new_face_1[2] = NODES[15];
  new_face_0[3] = NODES[13];  new_face_1[3] = NODES[17];

  new_face_2[0] = NODES[12];  new_face_3[0] = NODES[9];
  new_face_2[1] = NODES[11];  new_face_3[1] = NODES[12];
  new_face_2[2] = NODES[16];  new_face_3[2] = NODES[17];
  new_face_2[3] = NODES[17];  new_face_3[3] = NODES[14];

  for (i = 0; i < 4; i++)
  	compute_face_center(M, M->nfaces+i);
  
  // assemble children
  E_Int ncells = M->ncells;
  E_Int *new_cell_0 = &M->NFACE[6*ncells];
  E_Int *new_cell_1 = &M->NFACE[6*(ncells+1)];
  E_Int *new_cell_2 = &M->NFACE[6*(ncells+2)];
  E_Int *new_cell_3 = &M->NFACE[6*(ncells+3)];
  new_cell_0[0] = BOT[0];   new_cell_1[0] = BOT[1];
  new_cell_0[1] = TOP[0];   new_cell_1[1] = TOP[1];
  new_cell_0[2] = LFT[0];   new_cell_1[2] = nfaces;
  new_cell_0[3] = nfaces;   new_cell_1[3] = RGT[0];
  new_cell_0[4] = FRO[1];   new_cell_1[4] = FRO[0];
  new_cell_0[5] = nfaces+2; new_cell_1[5] = nfaces+3;
  
  new_cell_2[0] = BOT[3];   new_cell_3[0] = BOT[2];
  new_cell_2[1] = TOP[3];   new_cell_3[1] = TOP[2];
  new_cell_2[2] = LFT[1];   new_cell_3[2] = nfaces+1;
  new_cell_2[3] = nfaces+1; new_cell_3[3] = RGT[1];
  new_cell_2[4] = nfaces+2; new_cell_3[4] = nfaces+3;
  new_cell_2[5] = BCK[1];   new_cell_3[5] = BCK[0];
  
  for (i = 0; i < 4; i++)
  	compute_cell_center(M, M->ncells+i);
  
  // update cell tree
  tree_insert_children(ct, cell, M->ncells, 4);
  
  update_external_own_nei(cell, M, ct, ft);
  
  // set E_Internal face M->owner and M->neighbour
  M->owner[nfaces] = ncells;
  M->neigh[nfaces] = ncells+1;
  M->owner[nfaces+1] = ncells+2;
  M->neigh[nfaces+1] = ncells+3;
  M->owner[nfaces+2] = ncells;
  M->neigh[nfaces+2] = ncells+2;
  M->owner[nfaces+3] = ncells+1;
  M->neigh[nfaces+3] = ncells+3;

  M->ncells += 4;
  M->nfaces += 4;
}

void cut_cell_xz(E_Int cell, mesh *M, tree *ct, tree *ft)
{
  E_Int i, face, *pn, i0, reorient, *children, local[9];
  E_Int NODES[18] = {-1};
  E_Int FACES[16];
  E_Int *BOT = FACES;
  E_Int *TOP = FACES + 2;
  E_Int *LFT = FACES + 4;
  E_Int *RGT = FACES + 6;
  E_Int *FRO = FACES + 8;
  E_Int *BCK = FACES + 12;
  E_Int *pf = &M->NFACE[6*cell];
  E_Int *faces = M->NGON;
  
  // BOT
  face = pf[0];
  pn = &faces[4*face];
  i0 = 0;
  reorient = get_reorient(face, cell, normalIn[0], M);
  cut_face_x(face, M, ft, i0, reorient);
  assert(get_nchildren(ft, face) == 2);
  children = get_children(ft, face);
  for (i = 0; i < 2; i++) BOT[i] = children[i];
  Q6_get_ordered_data(face, M, i0, reorient, BOT, local);
  for (i = 0; i < 4; i++) NODES[i] = local[i];
  NODES[8] = local[4];
  NODES[13] = local[5];
  
  // LFT
  face = pf[2];
  pn = &faces[4*face];
  i0 = get_pos(NODES[0], pn, 4);
  reorient = get_reorient(face, cell, normalIn[2], M);
  cut_face_y(face, M, ft, i0, reorient);
  assert(get_nchildren(ft, face) == 2);
  children = get_children(ft, face);
  for (i = 0; i < 2; i++) LFT[i] = children[i];
  Q6_get_ordered_data(face, M, i0, reorient, LFT, local);
  assert(local[0] == NODES[0]);
  assert(local[1] == NODES[3]);
  NODES[4] = local[3];
  NODES[7] = local[2];
  NODES[16] = local[4];
  NODES[11] = local[5];
  
  // RGT
  face = pf[3];
  pn = &faces[4*face];
  i0 = get_pos(NODES[1], pn, 4);
  reorient = get_reorient(face, cell, normalIn[3], M);
  cut_face_y(face, M, ft, i0, reorient);
  assert(get_nchildren(ft, face) == 2);
  children = get_children(ft, face);
  for (i = 0; i < 2; i++) RGT[i] = children[i];
  Q6_get_ordered_data(face, M, i0, reorient, RGT, local);
  assert(local[0] == NODES[1]);
  assert(local[1] == NODES[2]);
  NODES[5] = local[3];
  NODES[6] = local[2];
  NODES[14] = local[4];
  NODES[9] = local[5];
  
  // TOP
  face = pf[1];
  pn = &faces[4*face];
  i0 = get_pos(NODES[4], pn, 4);
  reorient = get_reorient(face, cell, normalIn[1], M);
  cut_face_x(face, M, ft, i0, reorient);
  assert(get_nchildren(ft, face) == 2);
  children = get_children(ft, face);
  for (i = 0; i < 2; i++) TOP[i] = children[i];
  Q6_get_ordered_data(face, M, i0, reorient, TOP, local);
  assert(local[0] == NODES[4]);
  assert(local[1] == NODES[5]);
  assert(local[2] == NODES[6]);
  assert(local[3] == NODES[7]);
  NODES[10] = local[4];
  NODES[15] = local[5];
  
  // FRO
  face = pf[4];
  pn = &faces[4*face];
  i0 = get_pos(NODES[1], pn, 4);
  reorient = get_reorient(face, cell, normalIn[4], M);
  cut_face_xy(face, M, ft);
  assert(get_nchildren(ft, face) == 4);
  children = get_children(ft, face);
  for (i = 0; i < 4; i++) FRO[i] = children[i];
  Q9_get_ordered_data(face, M, i0, reorient, FRO, local);
  assert(local[0] == NODES[1]);
  assert(local[1] == NODES[0]);
  assert(local[2] == NODES[4]);
  assert(local[3] == NODES[5]);
  assert(local[4] == NODES[8]);
  assert(local[5] == NODES[11]);
  assert(local[6] == NODES[10]);
  assert(local[7] == NODES[9]);
  NODES[12] = local[8];
  
  // BCK
  face = pf[5];
  pn = &faces[4*face];
  i0 = get_pos(NODES[2], pn, 4);
  reorient = get_reorient(face, cell, normalIn[5], M);
  cut_face_xy(face, M, ft);
  assert(get_nchildren(ft, face) == 4);
  children = get_children(ft, face);
  for (i = 0; i < 4; i++) BCK[i] = children[i];
  Q9_get_ordered_data(face, M, i0, reorient, BCK, local);
  assert(local[0] == NODES[2]);
  assert(local[1] == NODES[3]);
  assert(local[2] == NODES[7]);
  assert(local[3] == NODES[6]);
  assert(local[4] == NODES[13]);
  assert(local[5] == NODES[16]);
  assert(local[6] == NODES[15]);
  assert(local[7] == NODES[14]);
  NODES[17] = local[8];
  
  // E_Internal faces
  E_Int nfaces = M->nfaces;
  E_Int *new_face_0 = &faces[4*nfaces];
  E_Int *new_face_1 = &faces[4*(nfaces+1)];
  E_Int *new_face_2 = &faces[4*(nfaces+2)];
  E_Int *new_face_3 = &faces[4*(nfaces+3)];
  new_face_0[0] = NODES[8];   new_face_1[0] = NODES[12];
  new_face_0[1] = NODES[13];  new_face_1[1] = NODES[17];
  new_face_0[2] = NODES[17];  new_face_1[2] = NODES[15];
  new_face_0[3] = NODES[12];  new_face_1[3] = NODES[10];

  new_face_2[0] = NODES[11];  new_face_3[0] = NODES[12];
  new_face_2[1] = NODES[12];  new_face_3[1] = NODES[9];
  new_face_2[2] = NODES[17];  new_face_3[2] = NODES[14];
  new_face_2[3] = NODES[16];  new_face_3[3] = NODES[17];

  for (i = 0; i < 4; i++)
  	compute_face_center(M, M->nfaces+i);
  
  // assemble children
  E_Int ncells = M->ncells;
  E_Int *new_cell_0 = &M->NFACE[6*ncells];
  E_Int *new_cell_1 = &M->NFACE[6*(ncells+1)];
  E_Int *new_cell_2 = &M->NFACE[6*(ncells+2)];
  E_Int *new_cell_3 = &M->NFACE[6*(ncells+3)];
  new_cell_0[0] = BOT[0];   new_cell_1[0] = BOT[1];
  new_cell_0[1] = nfaces+2; new_cell_1[1] = nfaces+3;
  new_cell_0[2] = LFT[0];   new_cell_1[2] = nfaces;
  new_cell_0[3] = nfaces;   new_cell_1[3] = RGT[0];
  new_cell_0[4] = FRO[1];   new_cell_1[4] = FRO[0];
  new_cell_0[5] = BCK[1];   new_cell_1[5] = BCK[0];
  new_cell_2[0] = nfaces+2; new_cell_3[0] = nfaces+3;
  new_cell_2[1] = TOP[0];   new_cell_3[1] = TOP[1];
  new_cell_2[2] = LFT[1];   new_cell_3[2] = nfaces+1;
  new_cell_2[3] = nfaces+1; new_cell_3[3] = RGT[1];
  new_cell_2[4] = FRO[2];   new_cell_3[4] = FRO[3];
  new_cell_2[5] = BCK[2];   new_cell_3[5] = BCK[3];
  
  for (i = 0; i < 4; i++)
  	compute_cell_center(M, M->ncells+i);
  
  // update cell tree
  tree_insert_children(ct, cell, M->ncells, 4);
  
  update_external_own_nei(cell, M, ct, ft);
  
  // set E_Internal face M->owner and M->neighbour
  M->owner[nfaces] = ncells;
  M->neigh[nfaces] = ncells+1;
  M->owner[nfaces+1] = ncells+2;
  M->neigh[nfaces+1] = ncells+3;
  M->owner[nfaces+2] = ncells;
  M->neigh[nfaces+2] = ncells+2;
  M->owner[nfaces+3] = ncells+1;
  M->neigh[nfaces+3] = ncells+3;
  
  M->ncells += 4;
  M->nfaces += 4;
}

void cut_cell_yz(E_Int cell, mesh *M, tree *ct, tree *ft)
{
  E_Int i, face, *pn, i0, reorient, *children, local[9];
  E_Int NODES[18] = {-1};
  E_Int FACES[16];
  E_Int *BOT = FACES;
  E_Int *TOP = FACES + 2;
  E_Int *LFT = FACES + 4;
  E_Int *RGT = FACES + 8;
  E_Int *FRO = FACES + 12;
  E_Int *BCK = FACES + 14;
  E_Int *pf = &M->NFACE[6*cell];
  E_Int *faces = M->NGON;
  
  // BOT
  face = pf[0];
  pn = &faces[4*face];
  i0 = 0;
  reorient = get_reorient(face, cell, normalIn[0], M);
  cut_face_y(face, M, ft, i0, reorient);
  assert(get_nchildren(ft, face) == 2);
  children = get_children(ft, face);
  for (i = 0; i < 2; i++) BOT[i] = children[i];
  Q6_get_ordered_data(face, M, i0, reorient, BOT, local);
  for (i = 0; i < 4; i++) NODES[i] = local[i];
  NODES[13] = local[4];
  NODES[8] = local[5];
  
  // LFT
  face = pf[2];
  pn = &faces[4*face];
  i0 = get_pos(NODES[0], pn, 4);
  reorient = get_reorient(face, cell, normalIn[2], M);
  cut_face_xy(face, M, ft);
  assert(get_nchildren(ft, face) == 4);
  children = get_children(ft, face);
  for (i = 0; i < 4; i++) LFT[i] = children[i];
  Q9_get_ordered_data(face, M, i0, reorient, LFT, local);
  assert(local[0] == NODES[0]);
  assert(local[1] == NODES[3]);
  assert(local[4] == NODES[8]);
  NODES[7] = local[2];
  NODES[4] = local[3];
  NODES[9] = local[5];
  NODES[10] = local[6];
  NODES[11] = local[7];
  NODES[12] = local[8];
  
  // RGT
  face = pf[3];
  pn = &faces[4*face];
  i0 = get_pos(NODES[1], pn, 4);
  reorient = get_reorient(face, cell, normalIn[3], M);
  cut_face_xy(face, M, ft);
  assert(get_nchildren(ft, face) == 4);
  children = get_children(ft, face);
  for (i = 0; i < 4; i++) RGT[i] = children[i];
  Q9_get_ordered_data(face, M, i0, reorient, RGT, local);
  assert(local[0] == NODES[1]);
  assert(local[1] == NODES[2]);
  assert(local[4] == NODES[13]);
  NODES[6] = local[2];
  NODES[5] = local[3];
  NODES[14] = local[5];
  NODES[15] = local[6];
  NODES[16] = local[7];
  NODES[17] = local[8];
  
  // TOP
  face = pf[1];
  pn = &faces[4*face];
  i0 = get_pos(NODES[4], pn, 4);
  reorient = get_reorient(face, cell, normalIn[1], M);
  cut_face_y(face, M, ft, i0, reorient);
  assert(get_nchildren(ft, face) == 2);
  children = get_children(ft, face);
  for (i = 0; i < 2; i++) TOP[i] = children[i];
  Q6_get_ordered_data(face, M, i0, reorient, TOP, local);
  assert(local[0] == NODES[4]);
  assert(local[1] == NODES[5]);
  assert(local[2] == NODES[6]);
  assert(local[3] == NODES[7]);
  assert(local[4] == NODES[15]);
  assert(local[5] == NODES[10]);
  
  // FRO
  face = pf[4];
  pn = &faces[4*face];
  i0 = get_pos(NODES[1], pn, 4);
  reorient = get_reorient(face, cell, normalIn[4], M);
  cut_face_y(face, M, ft, i0, reorient);
  assert(get_nchildren(ft, face) == 2);
  children = get_children(ft, face);
  for (i = 0; i < 2; i++) FRO[i] = children[i];
  Q6_get_ordered_data(face, M, i0, reorient, FRO, local);
  assert(local[0] == NODES[1]);
  assert(local[1] == NODES[0]);
  assert(local[2] == NODES[4]);
  assert(local[3] == NODES[5]);
  assert(local[4] == NODES[11]);
  assert(local[5] == NODES[16]);
  
  // BCK
  face = pf[5];
  pn = &faces[4*face];
  i0 = get_pos(NODES[2], pn, 4);
  reorient = get_reorient(face, cell, normalIn[5], M);
  cut_face_y(face, M, ft, i0, reorient);
  assert(get_nchildren(ft, face) == 2);
  children = get_children(ft, face);
  for (i = 0; i < 2; i++) BCK[i] = children[i];
  Q6_get_ordered_data(face, M, i0, reorient, BCK, local);
  assert(local[0] == NODES[2]);
  assert(local[1] == NODES[3]);
  assert(local[2] == NODES[7]);
  assert(local[3] == NODES[6]);
  assert(local[4] == NODES[9]);
  assert(local[5] == NODES[14]);
  
  // E_Internal faces
  E_Int nfaces = M->nfaces;
  E_Int *new_face_0 = &faces[4*nfaces];
  E_Int *new_face_1 = &faces[4*(nfaces+1)];
  E_Int *new_face_2 = &faces[4*(nfaces+2)];
  E_Int *new_face_3 = &faces[4*(nfaces+3)];
  new_face_0[0] = NODES[13];  new_face_1[0] = NODES[17];
  new_face_0[1] = NODES[8];   new_face_1[1] = NODES[12];
  new_face_0[2] = NODES[12];  new_face_1[2] = NODES[10];
  new_face_0[3] = NODES[17];  new_face_1[3] = NODES[15];

  new_face_2[0] = NODES[11];  new_face_3[0] = NODES[12];
  new_face_2[1] = NODES[16];  new_face_3[1] = NODES[17];
  new_face_2[2] = NODES[17];  new_face_3[2] = NODES[14];
  new_face_2[3] = NODES[12];  new_face_3[3] = NODES[9];

  for (i = 0; i < 4; i++)
  	compute_face_center(M, M->nfaces+i);
  
  // assemble children
  E_Int ncells = M->ncells;
  E_Int *new_cell_0 = &M->NFACE[6*ncells];
  E_Int *new_cell_1 = &M->NFACE[6*(ncells+1)];
  E_Int *new_cell_2 = &M->NFACE[6*(ncells+2)];
  E_Int *new_cell_3 = &M->NFACE[6*(ncells+3)];
  new_cell_0[0] = BOT[0];   new_cell_1[0] = BOT[1];
  new_cell_0[1] = nfaces+2; new_cell_1[1] = nfaces+3;
  new_cell_0[2] = LFT[0];   new_cell_1[2] = LFT[1];
  new_cell_0[3] = RGT[0];   new_cell_1[3] = RGT[1];
  new_cell_0[4] = FRO[0];   new_cell_1[4] = nfaces;
  new_cell_0[5] = nfaces;   new_cell_1[5] = BCK[0];
  
  new_cell_2[0] = nfaces+2; new_cell_3[0] = nfaces+3;
  new_cell_2[1] = TOP[0];   new_cell_3[1] = TOP[1];
  new_cell_2[2] = LFT[3];   new_cell_3[2] = LFT[2];
  new_cell_2[3] = RGT[3];   new_cell_3[3] = RGT[2];
  new_cell_2[4] = FRO[1];   new_cell_3[4] = nfaces+1;
  new_cell_2[5] = nfaces+1; new_cell_3[5] = BCK[1];
  for (i = 0; i < 4; i++)
  	compute_cell_center(M, M->ncells+i);
  
  // update cell tree
  tree_insert_children(ct, cell, M->ncells, 4);

  update_external_own_nei(cell, M, ct, ft);

  // set E_Internal face M->owner and M->neighbour
  M->owner[nfaces] = ncells;
  M->neigh[nfaces] = ncells+1;
  M->owner[nfaces+1] = ncells+2;
  M->neigh[nfaces+1] = ncells+3;
  M->owner[nfaces+2] = ncells;
  M->neigh[nfaces+2] = ncells+2;
  M->owner[nfaces+3] = ncells+1;
  M->neigh[nfaces+3] = ncells+3;

  M->ncells += 4;
  M->nfaces += 4;
}

void cut_cell_xyz(E_Int cell, mesh *M, tree *ct, tree *ft)
{
  E_Int i, face, *pn, i0, reorient, *children, local[9];
  E_Int NODES[27];
  E_Int FACES[24];
  E_Int *BOT = FACES;
  E_Int *TOP = FACES + 4;
  E_Int *LFT = FACES + 8;
  E_Int *RGT = FACES + 12;
  E_Int *FRO = FACES + 16;
  E_Int *BCK = FACES + 20;
  E_Int *pf = &M->NFACE[6*cell];
  E_Int *faces = M->NGON;
  
  // BOT
  face = pf[0];
  pn = &faces[4*face];
  i0 = 0;
  reorient = get_reorient(face, cell, normalIn[0], M);
  cut_face_xy(face, M, ft);
  assert(get_nchildren(ft, face) == 4);
  children = get_children(ft, face);
  for (i = 0; i < 4; i++) BOT[i] = children[i];
  Q9_get_ordered_data(face, M, i0, reorient, BOT, local);
  for (i = 0; i < 4; i++) NODES[i] = local[i];
  NODES[8] = local[4];
  NODES[9] = local[5];
  NODES[10] = local[6];
  NODES[11] = local[7];
  NODES[12] = local[8];
  
  // LFT
  face = pf[2];
  pn = &faces[4*face];
  i0 = get_pos(NODES[0], pn, 4);
  reorient = get_reorient(face, cell, normalIn[2], M);
  cut_face_xy(face, M, ft);
  assert(get_nchildren(ft, face) == 4);
  children = get_children(ft, face);
  for (i = 0; i < 4; i++) LFT[i] = children[i];
  Q9_get_ordered_data(face, M, i0, reorient, LFT, local);
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
  pn = &faces[4*face];
  i0 = get_pos(NODES[1], pn, 4);
  reorient = get_reorient(face, cell, normalIn[3], M);
  cut_face_xy(face, M, ft);
  assert(get_nchildren(ft, face) == 4);
  children = get_children(ft, face);
  for (i = 0; i < 4; i++) RGT[i] = children[i];
  Q9_get_ordered_data(face, M, i0, reorient, RGT, local);
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
  pn = &faces[4*face];
  i0 = get_pos(NODES[1], pn, 4);
  reorient = get_reorient(face, cell, normalIn[4], M);
  cut_face_xy(face, M, ft);
  assert(get_nchildren(ft, face) == 4);
  children = get_children(ft, face);
  for (i = 0; i < 4; i++) FRO[i] = children[i];
  Q9_get_ordered_data(face, M, i0, reorient, FRO, local);
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
  pn = &faces[4*face];
  i0 = get_pos(NODES[2], pn, 4);
  reorient = get_reorient(face, cell, normalIn[5], M);
  cut_face_xy(face, M, ft);
  assert(get_nchildren(ft, face) == 4);
  children = get_children(ft, face);
  for (i = 0; i < 4; i++) BCK[i] = children[i];
  Q9_get_ordered_data(face, M, i0, reorient, BCK, local);
  assert(local[0] == NODES[2]);
  assert(local[1] == NODES[3]);
  assert(local[2] == NODES[7]);
  assert(local[3] == NODES[6]);
  assert(local[4] == NODES[10]);
  assert(local[5] == NODES[13]);
  assert(local[7] == NODES[17]);
  NODES[23] = local[6];
  NODES[25] = local[8];
  
  // TOP
  face = pf[1];
  pn = &faces[4*face];
  i0 = get_pos(NODES[4], pn, 4);
  reorient = get_reorient(face, cell, normalIn[1], M);
  cut_face_xy(face, M, ft);
  assert(get_nchildren(ft, face) == 4);
  children = get_children(ft, face);
  for (i = 0; i < 4; i++) TOP[i] = children[i];
  Q9_get_ordered_data(face, M, i0, reorient, TOP, local);
  assert(local[0] == NODES[4]);
  assert(local[1] == NODES[5]);
  assert(local[2] == NODES[6]);
  assert(local[3] == NODES[7]);
  assert(local[4] == NODES[21]);
  assert(local[5] == NODES[18]);
  assert(local[6] == NODES[23]);
  assert(local[7] == NODES[14]);
  NODES[24] = local[8];
  
  // Add cell centroid
  NODES[26] = M->npoints;
  compute_cell_center(M, cell);
  E_Float *cx = &M->cc[3*cell];
  E_Float *X = &M->xyz[3*M->npoints];
  for (E_Int k = 0; k < 3; k++)
  	X[k] = cx[k];
  M->npoints++;
  
  // E_Internal faces
  E_Int nfaces = M->nfaces;
  E_Int *new_face_0 = &faces[4*nfaces];
  E_Int *new_face_1 = &faces[4*(nfaces+1)];
  E_Int *new_face_2 = &faces[4*(nfaces+2)];
  E_Int *new_face_3 = &faces[4*(nfaces+3)];
  E_Int *new_face_4 = &faces[4*(nfaces+4)];
  E_Int *new_face_5 = &faces[4*(nfaces+5)];
  E_Int *new_face_6 = &faces[4*(nfaces+6)];
  E_Int *new_face_7 = &faces[4*(nfaces+7)];
  E_Int *new_face_8 = &faces[4*(nfaces+8)];
  E_Int *new_face_9 = &faces[4*(nfaces+9)];
  E_Int *new_face_10 = &faces[4*(nfaces+10)];
  E_Int *new_face_11 = &faces[4*(nfaces+11)];
  new_face_0[0] = NODES[8];   new_face_1[0] = NODES[12];
  new_face_0[1] = NODES[12];  new_face_1[1] = NODES[10];
  new_face_0[2] = NODES[26];  new_face_1[2] = NODES[25];
  new_face_0[3] = NODES[22];  new_face_1[3] = NODES[26];
  
  new_face_2[0] = NODES[15];  new_face_3[0] = NODES[22];
  new_face_2[1] = NODES[22];  new_face_3[1] = NODES[19];
  new_face_2[2] = NODES[26];  new_face_3[2] = NODES[20];
  new_face_2[3] = NODES[16];  new_face_3[3] = NODES[26];

  new_face_4[0] = NODES[16];  new_face_5[0] = NODES[26];
  new_face_4[1] = NODES[26];  new_face_5[1] = NODES[20];
  new_face_4[2] = NODES[25];  new_face_5[2] = NODES[17];
  new_face_4[3] = NODES[13];  new_face_5[3] = NODES[25];
  
  new_face_6[0] = NODES[12];  new_face_7[0] = NODES[9]; 
  new_face_6[1] = NODES[11];  new_face_7[1] = NODES[12];
  new_face_6[2] = NODES[16];  new_face_7[2] = NODES[26];
  new_face_6[3] = NODES[26];  new_face_7[3] = NODES[20];
  
  new_face_8[0] = NODES[22];  new_face_9[0] = NODES[26];
  new_face_8[1] = NODES[26];  new_face_9[1] = NODES[25];
  new_face_8[2] = NODES[24];  new_face_9[2] = NODES[23];
  new_face_8[3] = NODES[21];  new_face_9[3] = NODES[24];

  new_face_10[0] = NODES[26];  new_face_11[0] = NODES[20];
  new_face_10[1] = NODES[16];  new_face_11[1] = NODES[26];
  new_face_10[2] = NODES[14];  new_face_11[2] = NODES[24];
  new_face_10[3] = NODES[24];  new_face_11[3] = NODES[18];
  
  for (i = 0; i < 12; i++)
  	compute_face_center(M, M->nfaces+i);

  // assemble children
  E_Int ncells = M->ncells;
  E_Int *new_cell_0 = &M->NFACE[6*ncells];
  E_Int *new_cell_1 = &M->NFACE[6*(ncells+1)];
  E_Int *new_cell_2 = &M->NFACE[6*(ncells+2)];
  E_Int *new_cell_3 = &M->NFACE[6*(ncells+3)];
  E_Int *new_cell_4 = &M->NFACE[6*(ncells+4)];
  E_Int *new_cell_5 = &M->NFACE[6*(ncells+5)];
  E_Int *new_cell_6 = &M->NFACE[6*(ncells+6)];
  E_Int *new_cell_7 = &M->NFACE[6*(ncells+7)];
  new_cell_0[0] = BOT[0];   new_cell_1[0] = BOT[1];
  new_cell_0[1] = nfaces+2; new_cell_1[1] = nfaces+3;
  new_cell_0[2] = LFT[0];   new_cell_1[2] = nfaces;
  new_cell_0[3] = nfaces;   new_cell_1[3] = RGT[0];
  new_cell_0[4] = FRO[1];   new_cell_1[4] = FRO[0];
  new_cell_0[5] = nfaces+6; new_cell_1[5] = nfaces+7;
  
  new_cell_2[0] = BOT[3];   new_cell_3[0] = BOT[2];
  new_cell_2[1] = nfaces+4; new_cell_3[1] = nfaces+5;
  new_cell_2[2] = LFT[1];   new_cell_3[2] = nfaces+1;
  new_cell_2[3] = nfaces+1; new_cell_3[3] = RGT[1];
  new_cell_2[4] = nfaces+6; new_cell_3[4] = nfaces+7;
  new_cell_2[5] = BCK[1];   new_cell_3[5] = BCK[0];
  
  new_cell_4[0] = nfaces+2;  new_cell_5[0] = nfaces+3;
  new_cell_4[1] = TOP[0];    new_cell_5[1] = TOP[1];
  new_cell_4[2] = LFT[3];    new_cell_5[2] = nfaces+8;
  new_cell_4[3] = nfaces+8;  new_cell_5[3] = RGT[3];
  new_cell_4[4] = FRO[2];    new_cell_5[4] = FRO[3];
  new_cell_4[5] = nfaces+10; new_cell_5[5] = nfaces+11;
  
  new_cell_6[0] = nfaces+4;  new_cell_7[0] = nfaces+5;
  new_cell_6[1] = TOP[3];    new_cell_7[1] = TOP[2];
  new_cell_6[2] = LFT[2];    new_cell_7[2] = nfaces+9;
  new_cell_6[3] = nfaces+9;  new_cell_7[3] = RGT[2];
  new_cell_6[4] = nfaces+10; new_cell_7[4] = nfaces+11;
  new_cell_6[5] = BCK[2];    new_cell_7[5] = BCK[3];

  //for (i = 0; i < 8; i++)
  //	compute_cell_center(M, M->ncells+i);
  
  // update cell tree
  tree_insert_children(ct, cell, M->ncells, 8);
  
  // set external faces M->owner and M->neighbours
  update_external_own_nei(cell, M, ct, ft);
  
  // set E_Internal faces M->owners and M->neighbours
  M->owner[nfaces] = ncells;      M->owner[nfaces+1] = ncells+2;
  M->neigh[nfaces] = ncells+1;    M->neigh[nfaces+1] = ncells+3;

  M->owner[nfaces+2] = ncells;    M->owner[nfaces+3] = ncells+1;
  M->neigh[nfaces+2] = ncells+4;  M->neigh[nfaces+3] = ncells+5;

  M->owner[nfaces+4] = ncells+2;  M->owner[nfaces+5] = ncells+3;
  M->neigh[nfaces+4] = ncells+6;  M->neigh[nfaces+5] = ncells+7;

  M->owner[nfaces+6] = ncells;    M->owner[nfaces+7] = ncells+1;
  M->neigh[nfaces+6] = ncells+2;  M->neigh[nfaces+7] = ncells+3;
  
  M->owner[nfaces+8] = ncells+4;  M->owner[nfaces+9] = ncells+6;
  M->neigh[nfaces+8] = ncells+5;  M->neigh[nfaces+9] = ncells+7;

  M->owner[nfaces+10] = ncells+4; M->owner[nfaces+11] = ncells+5;
  M->neigh[nfaces+10] = ncells+6; M->neigh[nfaces+11] = ncells+7;

  M->ncells += 8;
  M->nfaces += 12;
}
