#include "Proto.h"

void Q9_unrefine(E_Int face, AMesh *M)
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
}

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

void Q9_refine(E_Int face, AMesh *M)
{
  E_Int *pn = get_facets(face, M->ngon, M->indPG);

  Edge E;
  E_Int ec[4]; // Edge center nodes

  for (E_Int i = 0; i < 4; i++) {
    E_Int ni = pn[i];
    E_Int nj = pn[(i+1)%4];
    E.set(ni, nj);

    auto search = M->ecenter->find(E);
    if (search == M->ecenter->end()) {
      M->x[M->npoints] = 0.5*(M->x[ni] + M->x[nj]);
      M->y[M->npoints] = 0.5*(M->y[ni] + M->y[nj]);
      M->z[M->npoints] = 0.5*(M->z[ni] + M->z[nj]);
      ec[i] = M->npoints;
      (*M->ecenter)[E] = M->npoints++;
    } else {
      ec[i] = search->second;
    }
  }

  // Note(Imad): face center is supposed to be already computed
  E_Float fc[3] = {0., 0., 0.};
  for (E_Int i = 0; i < 4; i++) {
    fc[0] += M->x[pn[i]];
    fc[1] += M->y[pn[i]];
    fc[2] += M->z[pn[i]];
  }
  for (E_Int i = 0; i < 3; i++) fc[i] /= 4.0;
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
  E_Int next_level = M->faceTree->level(face) + 1;
  M->faceTree->set_parent_elem(face, 4, M->nfaces);
  for (E_Int i = 1; i < 4; i++)
    M->faceTree->set_child_elem(i, face, QUAD, next_level, M->nfaces);

  M->owner[M->nfaces] = M->owner[face];
  M->neigh[M->nfaces] = M->neigh[face];
  M->owner[M->nfaces+1] = M->owner[face];
  M->neigh[M->nfaces+1] = M->neigh[face];
  M->owner[M->nfaces+2] = M->owner[face];
  M->neigh[M->nfaces+2] = M->neigh[face];

  M->nfaces += 3;
}

void Q6_refine(E_Int face, AMesh *M)
{
  E_Int *pn = get_facets(face, M->ngon, M->indPG);

  Edge E;
  E_Int ec[2];

  E_Int ni, nj;

  // First first edge that's perpendicular to 2D normal vector

  E_Int start = 0;
  for (; start < 4; start++) {
    ni = pn[start];
    nj = pn[(start+1)%4];

    E_Float e[3] = {M->x[nj]-M->x[ni],
                    M->y[nj]-M->y[ni],
                    M->z[nj]-M->z[ni]};

    E_Float dp = K_MATH::dot(e, M->mode_2D, 3);

    if (fabs(dp) < 1.0e-13) break;
  }

  assert(start < 4);

  // First edge
  ni = pn[start];
  nj = pn[(start+1)%4];
  E.set(ni, nj);

  {
    auto search = M->ecenter->find(E);
    if (search == M->ecenter->end()) {
      M->x[M->npoints] = 0.5*(M->x[ni] + M->x[nj]);
      M->y[M->npoints] = 0.5*(M->y[ni] + M->y[nj]);
      M->z[M->npoints] = 0.5*(M->z[ni] + M->z[nj]);
      ec[0] = M->npoints;
      (*M->ecenter)[E] = M->npoints++;
    } else {
      ec[0] = search->second;
    }
  }

  // Second edge
  ni = pn[(start+2)%4];
  nj = pn[(start+3)%4];
  E.set(ni, nj);

  {
    auto search = M->ecenter->find(E);
    if (search == M->ecenter->end()) {
      M->x[M->npoints] = 0.5*(M->x[ni] + M->x[nj]);
      M->y[M->npoints] = 0.5*(M->y[ni] + M->y[nj]);
      M->z[M->npoints] = 0.5*(M->z[ni] + M->z[nj]);
      ec[1] = M->npoints;
      (*M->ecenter)[E] = M->npoints++;
    } else {
      ec[1] = search->second;
    }
  }

  // Set the rest of the children in ngon, starting from ptr
  E_Int *ptr = &M->indPG[M->nfaces];

  if (start == 0 || start == 2) {
    // dir X
    
    // Second child
    ptr[1] = *ptr + 4;
    M->ngon[*ptr    ] = ec[0];
    M->ngon[*ptr + 1] = pn[1];
    M->ngon[*ptr + 2] = pn[2];
    M->ngon[*ptr + 3] = ec[1];

    // First child replaces face
    ptr = &M->indPG[face];
    //M->ngon[*ptr    ] = pn[start];
    M->ngon[*ptr + 1] = ec[0];
    M->ngon[*ptr + 2] = ec[1];
    //M->ngon[*ptr + 3] = ec[1];
  } else {
    
    // Second child
    ptr[1] = *ptr + 4;
    M->ngon[*ptr    ] = ec[1];
    M->ngon[*ptr + 1] = ec[0];
    M->ngon[*ptr + 2] = pn[2];
    M->ngon[*ptr + 3] = pn[3];

    // First child replaces face
    ptr = &M->indPG[face];
    //M->ngon[*ptr    ] = pn[start];
    //M->ngon[*ptr + 1] = ec[0];
    M->ngon[*ptr + 2] = ec[0];
    M->ngon[*ptr + 3] = ec[1];
  }

  // Set faces in faceTree
  E_Int next_level = M->faceTree->level(face) + 1;
  M->faceTree->set_parent_elem(face, 2, M->nfaces);
  for (E_Int i = 1; i < 2; i++)
    M->faceTree->set_child_elem(i, face, QUAD, next_level, M->nfaces);

  M->owner[M->nfaces] = M->owner[face];
  M->neigh[M->nfaces] = M->neigh[face];

  M->nfaces += 1;
}

void Q6_unrefine(E_Int face, AMesh *M)
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

  E_Int dir = pn0[1] == pn1[0] ? DIRX : DIRY; 

  if (dir == DIRX) {
    pn0[1] = pn1[1];
    pn0[2] = pn1[2];
  } else {
    pn0[2] = pn1[2];
    pn0[3] = pn1[3];
  }
}


void Q6_get_ordered_data(AMesh *M, E_Int NODE, E_Int reorient,
  E_Int *children, E_Int *local)
{
  E_Int *pn0 = get_facets(children[0], M->ngon, M->indPG);
  E_Int *pn1 = get_facets(children[1], M->ngon, M->indPG);

  E_Int dir = pn0[1] == pn1[0] ? DIRX : DIRY;

  if (dir == DIRX) {
    local[0] = pn0[0];
    local[1] = pn1[1];
    local[2] = pn1[2];
    local[3] = pn0[3];

    E_Int i0 = Get_pos(NODE, local, 4);
    assert(i0 != -1);
  
    Right_shift(local, i0, 4);

    if (!reorient) {
      switch (i0) {
        case 0:
        case 3:
          local[4] = pn0[1];
          local[5] = pn0[2];
          break;
        case 1:
        case 2:
          local[4] = pn0[2];
          local[5] = pn0[1];
          std::swap(children[0], children[1]);
          break;
        default:
          assert(0);
          break;
      }
    } else {
      std::swap(local[1], local[3]);
      
      switch (i0) {
        case 0:
        case 3:
          local[4] = pn0[2];
          local[5] = pn0[1];
          break;
        case 1:
        case 2:
          local[4] = pn0[1];
          local[5] = pn0[2];
          std::swap(children[0], children[1]);
          break;
        default:
          assert(0);
          break;
      }
    }
  } else { // dir = dirY
    
    local[0] = pn0[0];
    local[1] = pn0[1];
    local[2] = pn1[2];
    local[3] = pn1[3];

    E_Int i0 = Get_pos(NODE, local, 4);
    assert(i0 != -1);
  
    Right_shift(local, i0, 4);

    if (!reorient) {
      switch (i0) {
        case 0:
        case 1:
          local[4] = pn0[2];
          local[5] = pn0[3];
          break;
        case 2:
        case 3:
          local[4] = pn0[3];
          local[5] = pn0[2];
          std::swap(children[0], children[1]);
          break;
        default:
          assert(0);
          break;
      }
    } else {
      std::swap(local[1], local[3]);
      
      switch (i0) {
        case 0:
        case 1:
          local[4] = pn0[3];
          local[5] = pn0[2];
          break;
        case 2:
        case 3:
          local[4] = pn0[2];
          local[5] = pn0[3];
          std::swap(children[0], children[1]);
          break;
        default:
          assert(0);
          break;
      }
    }
  }
}
