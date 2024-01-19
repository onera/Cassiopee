#include "Proto.h"

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

void refine_quad(E_Int face, AMesh *M)
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
