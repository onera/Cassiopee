#include "Proto.h"

void get_ref_cells_and_faces(AMesh *M, std::vector<E_Int> &ref_cells,
  std::vector<E_Int> &ref_faces)
{
  ref_cells.clear();
  ref_faces.clear();

  E_Int *vcells = (E_Int *)XCALLOC(M->ncells, sizeof(E_Int));

  for (E_Int i = 0; i < M->nfaces; i++) {
    E_Int own = M->owner[i];
    if (M->ref_data[3*own] > 0) {
      ref_faces.push_back(i);
      if (!vcells[own]) {
        vcells[own] = 1;
        ref_cells.push_back(own);
      }
      continue;
    }

    E_Int nei = M->neigh[i];
    if (nei == -1) continue;

    if (M->ref_data[3*nei] > 0) {
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
void set_face(Element **tree, E_Int where, E_Int id, E_Int nchildren,
  E_Int parent, E_Int position, E_Int type, E_Int level)
{
  if (nchildren > 0)
    tree[id]->children = (E_Int *)XMALLOC(nchildren * sizeof(E_Int));
  
  for (E_Int i = 0; i < nchildren; i++)
    tree[id]->children[i] = where + i;

  tree[id]->nchildren = nchildren;
  tree[id]->parent    = parent;
  tree[id]->position  = position;
  tree[id]->type      = type;
  tree[id]->level     = level;
}

static
void set_cell(Element **tree, E_Int where, E_Int id, E_Int nchildren,
  E_Int parent, E_Int position, E_Int type, E_Int level)
{
  if (nchildren > 0)
    tree[id]->children = (E_Int *)XMALLOC(nchildren * sizeof(E_Int));
  
  for (E_Int i = 0; i < nchildren; i++)
    tree[id]->children[i] = where + i;

  tree[id]->nchildren = nchildren;
  tree[id]->parent    = parent;
  tree[id]->position  = position;
  tree[id]->type      = type;
  tree[id]->level     = level;
}

static
void refine_tri(E_Int face, AMesh *M)
{
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

  E_Int next_level = M->faceTree[face]->level + 1;

  // First child replaces parent in faceTree and ngon

  // Set faces in faceTree
  set_face(M->faceTree, M->nfaces, face,        3, face, 0, TRI, next_level);
  set_face(M->faceTree, M->nfaces, M->nfaces,   0, face, 1, TRI, next_level);
  set_face(M->faceTree, M->nfaces, M->nfaces+1, 0, face, 2, TRI, next_level);
  set_face(M->faceTree, M->nfaces, M->nfaces+2, 0, face, 3, TRI, next_level);

  // Set faces in ngon, starting from ptr
  E_Int *ptr = &M->indPG[M->nfaces];

  // Second face
  ptr[1] = *ptr + 3;
  M->ngon[*ptr]     = pn[1];
  M->ngon[*ptr + 1] = ec[1];
  M->ngon[*ptr + 2] = ec[0];
  ptr++;

  // Third face
  ptr[1] = *ptr + 3;
  M->ngon[*ptr]     = pn[2];
  M->ngon[*ptr + 1] = ec[2];
  M->ngon[*ptr + 2] = ec[1];
  ptr++;

  // Fourth face
  ptr[1] = *ptr + 3;
  M->ngon[*ptr] = pn[2];
  M->ngon[*ptr + 1] = ec[2];
  M->ngon[*ptr + 2] = ec[1];

  // First face
  pn[1] = ec[0];
  pn[2] = ec[2];

  M->owner[M->nfaces] = M->owner[face];
  M->neigh[M->nfaces] = M->neigh[face];
  M->owner[M->nfaces+1] = M->owner[face];
  M->neigh[M->nfaces+1] = M->neigh[face];
  M->owner[M->nfaces+2] = M->owner[face];
  M->neigh[M->nfaces+2] = M->neigh[face];

  M->nfaces += 3;
}

static
void refine_quad(E_Int face, AMesh *M)
{
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

  E_Float fc[3], fa[3];
  compute_face_center_and_area(face, 4, pn, M->x, M->y, M->z, fc, fa);
  M->x[M->npoints] = fc[0];
  M->y[M->npoints] = fc[1];
  M->z[M->npoints] = fc[2];
  E_Int ncenter = M->npoints++;

  E_Int next_level = M->faceTree[face]->level + 1;

  // First child replaces parent in faceTree and ngon

  // Set faces in faceTree
  set_face(M->faceTree, M->nfaces, face,        3, face, 0, QUAD, next_level);
  set_face(M->faceTree, M->nfaces, M->nfaces,   0, face, 1, QUAD, next_level);
  set_face(M->faceTree, M->nfaces, M->nfaces+1, 0, face, 2, QUAD, next_level);
  set_face(M->faceTree, M->nfaces, M->nfaces+2, 0, face, 3, QUAD, next_level);

  // Set faces in ngon, starting from ptr
  E_Int *ptr = &M->indPG[M->nfaces];

  // Second face
  ptr[1] = *ptr + 4;
  M->ngon[*ptr]     = pn[1];
  M->ngon[*ptr + 1] = ec[1];
  M->ngon[*ptr + 2] = ncenter;
  M->ngon[*ptr + 3] = ec[0];
  ptr++;

  // Third face
  ptr[1] = *ptr + 4;
  M->ngon[*ptr]     = pn[2];
  M->ngon[*ptr + 1] = ec[2];
  M->ngon[*ptr + 2] = ncenter;
  M->ngon[*ptr + 3] = ec[1];
  ptr++;

  // Fourth face
  ptr[1] = *ptr + 4;
  M->ngon[*ptr] = pn[3];
  M->ngon[*ptr + 1] = ec[3];
  M->ngon[*ptr + 2] = ncenter;
  M->ngon[*ptr + 3] = ec[2];

  // First face
  pn[1] = ec[0];
  pn[2] = ncenter;
  pn[3] = ec[3];

  M->owner[M->nfaces] = M->owner[face];
  M->neigh[M->nfaces] = M->neigh[face];
  M->owner[M->nfaces+1] = M->owner[face];
  M->neigh[M->nfaces+1] = M->neigh[face];
  M->owner[M->nfaces+2] = M->owner[face];
  M->neigh[M->nfaces+2] = M->neigh[face];

  M->nfaces += 3;
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

static
void Q9_get_ordered_data(E_Int i0, E_Int reorient, E_Int *children,
  E_Int *pn)
{
  Right_shift(pn, i0, 4);
  Right_shift(pn+4, i0, 4);
  Right_shift(children, i0, 4);

  if (reorient) {
  	std::swap(pn[1], pn[3]);
  	std::swap(pn[4], pn[7]);
  	std::swap(pn[5], pn[6]);
  	std::swap(children[1], children[3]);
  }
}

static
void reconstruct_parent_face(E_Int face, AMesh *M, E_Int *children, E_Int *pn)
{
  E_Int *ptr = M->faceTree[face]->children;
  children[0] = face;
  children[1] = ptr[0];
  children[2] = ptr[1];
  children[3] = ptr[2];

  E_Int *pn0 = get_facets(children[0], M->ngon, M->indPG);
  E_Int *pn1 = get_facets(children[1], M->ngon, M->indPG);
  E_Int *pn2 = get_facets(children[2], M->ngon, M->indPG);
  E_Int *pn3 = get_facets(children[3], M->ngon, M->indPG);

  pn[0] = pn0[0];
  pn[1] = pn1[0];
  pn[2] = pn2[0];
  pn[3] = pn3[0];
  
  pn[4] = pn0[1];
  pn[5] = pn1[1];
  pn[6] = pn2[1];
  pn[7] = pn3[1];
 
  pn[8] = pn0[2];
}

static
void refine_hexa(E_Int cell, AMesh *M)
{
  E_Int NODES[27], FACES[24];
  E_Int *BOT = FACES;
  E_Int *TOP = FACES + 4;
  E_Int *LFT = FACES + 8;
  E_Int *RGT = FACES + 12;
  E_Int *FRO = FACES + 16;
  E_Int *BCK = FACES + 20;
  E_Int *pf = get_facets(cell, M->nface, M->indPH);

  E_Int face, i0, reorient, pn[9];

  // BOT
  face = pf[0];
  reconstruct_parent_face(face, M, BOT, pn);
  i0 = 0;
  reorient = get_reorient(face, cell, normalIn_H[0], M);
  Q9_get_ordered_data(i0, reorient, BOT, pn);
  for (E_Int i = 0; i < 4; i++) NODES[i] = pn[i];
  NODES[8] = pn[4];
  NODES[9] = pn[5];
  NODES[10] = pn[6];
  NODES[11] = pn[7];
  NODES[12] = pn[8];

  // LFT
  face = pf[2];
  reconstruct_parent_face(face, M, LFT, pn);
  i0 = Get_pos(NODES[0], pn, 4);
  reorient = get_reorient(face, cell, normalIn_H[2], M);
  Q9_get_ordered_data(i0, reorient, LFT, pn);
  assert(pn[0] == NODES[0]);
  assert(pn[1] == NODES[3]);
  assert(pn[4] == NODES[11]);
  NODES[7] = pn[2];
  NODES[4] = pn[3];
  NODES[13] = pn[5];
  NODES[14] = pn[6];
  NODES[15] = pn[7];
  NODES[16] = pn[8];

  // RGT
  face = pf[3];
  reconstruct_parent_face(face, M, RGT, pn);
  i0 = Get_pos(NODES[1], pn, 4);
  reorient = get_reorient(face, cell, normalIn_H[3], M);
  Q9_get_ordered_data(i0, reorient, RGT, pn);
  assert(pn[0] == NODES[1]);
  assert(pn[1] == NODES[2]);
  assert(pn[4] == NODES[9]);
  NODES[6] = pn[2];
  NODES[5] = pn[3];
  NODES[17] = pn[5];
  NODES[18] = pn[6];
  NODES[19] = pn[7];
  NODES[20] = pn[8];

  // FRO
  face = pf[4];
  reconstruct_parent_face(face, M, FRO, pn);
  i0 = Get_pos(NODES[1], pn, 4);
  reorient = get_reorient(face, cell, normalIn_H[4], M);
  Q9_get_ordered_data(i0, reorient, FRO, pn);
  assert(pn[0] == NODES[1]);
  assert(pn[1] == NODES[0]);
  assert(pn[2] == NODES[4]);
  assert(pn[3] == NODES[5]);
  assert(pn[4] == NODES[8]);
  assert(pn[5] == NODES[15]);
  assert(pn[7] == NODES[19]);
  NODES[21] = pn[6];
  NODES[22] = pn[8];

  // BCK
  face = pf[5];
  reconstruct_parent_face(face, M, BCK, pn);
  i0 = Get_pos(NODES[2], pn, 4);
  reorient = get_reorient(face, cell, normalIn_H[5], M);
  Q9_get_ordered_data(i0, reorient, BCK, pn);
  assert(pn[0] == NODES[2]);
  assert(pn[1] == NODES[3]);
  assert(pn[2] == NODES[7]);
  assert(pn[3] == NODES[6]);
  assert(pn[4] == NODES[10]);
  assert(pn[5] == NODES[13]);
  assert(pn[7] == NODES[17]);
  NODES[23] = pn[6];
  NODES[25] = pn[8];

  // TOP
  face = pf[1];
  reconstruct_parent_face(face, M, TOP, pn);
  i0 = Get_pos(NODES[4], pn, 4);
  reorient = get_reorient(face, cell, normalIn_H[1], M);
  Q9_get_ordered_data(i0, reorient, TOP, pn);
  assert(pn[0] == NODES[4]);
  assert(pn[1] == NODES[5]);
  assert(pn[2] == NODES[6]);
  assert(pn[3] == NODES[7]);
  assert(pn[4] == NODES[21]);
  assert(pn[5] == NODES[18]);
  assert(pn[6] == NODES[23]);
  assert(pn[7] == NODES[14]);
  NODES[24] = pn[8];
}

static
void refine_tetra(E_Int cell, AMesh *M)
{}

static
void refine_penta(E_Int cell, AMesh *M)
{}

static
void refine_pyra(E_Int cell, AMesh *M)
{}

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

void resize_data_for_refinement(AMesh *M, E_Int nref_cells, E_Int nref_faces)
{
  E_Int cell_incr = nref_cells * 7;
  E_Int face_incr = nref_faces * 3; // OK for quads and tris
  E_Int nnew_cells = M->ncells + cell_incr;
  E_Int nnew_faces = M->nfaces + face_incr;
  // max estimation: 5 new points per refined quad + nref_cells centroids
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
}
