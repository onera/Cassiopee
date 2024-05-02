#include "Proto.h"

void T6_refine(E_Int face, AMesh *M)
{
  E_Int *pn = get_facets(face, M->ngon, M->indPG);

  Edge E;
  E_Int ec[3]; // Edge center nodes

  for (E_Int i = 0; i < 3; i++) {
    E_Int ni = pn[i];
    E_Int nj = pn[(i+1)%3];
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
  E_Int next_level = M->faceTree->level(face) + 1;
  M->faceTree->set_parent_elem(face, 4, M->nfaces);
  for (E_Int i = 1; i < 4; i++)
    M->faceTree->set_child_elem(i, face, TRI, next_level, M->nfaces);

  // Own and nei
  M->owner[M->nfaces] = M->owner[face];
  M->neigh[M->nfaces] = M->neigh[face];
  M->owner[M->nfaces+1] = M->owner[face];
  M->neigh[M->nfaces+1] = M->neigh[face];
  M->owner[M->nfaces+2] = M->owner[face];
  M->neigh[M->nfaces+2] = M->neigh[face];

  M->nfaces += 3;
}

void T6_unrefine(E_Int face, AMesh *)
{}

void reorder_tri(E_Int *local, E_Int *pn, E_Int reorient, E_Int i0)
{
  for (E_Int i = 0; i < 3; i++) local[i] = pn[i];
  Right_shift(local, i0, 3);
  if (reorient) std::swap(local[1], local[2]);
}

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