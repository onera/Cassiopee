#include "Proto.h"

#define MINREF -10
#define RAISE(error) PyErr_SetString(PyExc_ValueError, (error))

PyObject *K_XCORE::CreateAdaptMesh(PyObject *self, PyObject *args)
{
  PyObject *ARRAY, *COMM;
  E_Float Tr;
  if (!PYPARSETUPLE_(args, OO_ R_, &ARRAY, &COMM, &Tr)) {
    RAISE("Wrong input.");
    return NULL;
  }

  // Check input
  E_Int ret;
  E_Int ni, nj, nk;
  K_FLD::FldArrayF *f;
  K_FLD::FldArrayI *cn;
  char *varString, *eltType;
  ret = K_ARRAY::getFromArray3(ARRAY, varString, f, ni, nj, nk, cn, eltType);

  if (ret <= 0) {
    RAISE("Bad mesh.");
    return NULL;
  }

  if (ret == 1) {
    RAISE("Mesh is not NGon.");
    RELEASESHAREDS(ARRAY, f);
    return NULL;
  }

  E_Int posx = K_ARRAY::isCoordinateXPresent(varString);
  E_Int posy = K_ARRAY::isCoordinateYPresent(varString);
  E_Int posz = K_ARRAY::isCoordinateZPresent(varString);

  if (posx == -1 || posy == -1 || posz == -1) {
    RELEASESHAREDU(ARRAY, f, cn);
    RAISE("Bad coordinates.");
    return NULL;
  }

  posx++; posy++; posz++;

  // Check mesh validity
  ret = K_CONNECT::check_overlapping_cells(*cn);
  if (ret != 0) {
    RELEASESHAREDU(ARRAY, f, cn);
    RAISE("Cells should not overlap.");
    return NULL;
  }

  ret = K_CONNECT::check_open_cells(*cn, NULL);
  if (ret != 0) {
    RELEASESHAREDU(ARRAY, f, cn);
    RAISE("Cells should be closed.");
    return NULL;
  }

  // Apply consistent orientation on mesh faces
  E_Float *px = f->begin(posx);
  E_Float *py = f->begin(posy);
  E_Float *pz = f->begin(posz);

  ret = K_CONNECT::orient_boundary_ngon(px, py, pz, *cn);

  if (ret != 0) {
    RELEASESHAREDU(ARRAY, f, cn);
    RAISE("Failed boundary orientation.");
    return NULL;
  }

  // Init mesh
  AMesh *M = new AMesh();
  M->nfaces = cn->getNFaces();
  M->owner = (E_Int *)XMALLOC(M->nfaces * sizeof(E_Int));
  M->neigh = (E_Int *)XMALLOC(M->nfaces * sizeof(E_Int));

  // Orient the internal faces and build parent elements
  K_CONNECT::build_parent_elements_ngon(*cn, M->owner, M->neigh);

  E_Int *ngon = cn->getNGon();
  E_Int *indPG = cn->getIndPG();
  E_Int *nface = cn->getNFace();
  E_Int *indPH = cn->getIndPH();
  M->ncells = cn->getNElts();
  M->npoints = f->getSize();

  M->indPH = (E_Int *)XMALLOC((M->ncells+1) * sizeof(E_Int));
  for (E_Int i = 0; i < M->ncells+1; i++) M->indPH[i] = indPH[i];

  M->indPG = (E_Int *)XMALLOC((M->nfaces+1) * sizeof(E_Int));
  for (E_Int i = 0; i < M->nfaces+1; i++) M->indPG[i] = indPG[i];

  M->nface = (E_Int *)XMALLOC(M->indPH[M->ncells] * sizeof(E_Int));
  for (E_Int i = 0; i < M->ncells; i++) {
    E_Int nf = -1;
    E_Int *pf = cn->getElt(i, nf, nface, indPH);
    E_Int *ptr = &M->nface[M->indPH[i]];
    for (E_Int j = 0; j < nf; j++)
      ptr[j] = pf[j]-1;
  }

  M->ngon = (E_Int *)XMALLOC(M->indPG[M->nfaces] * sizeof(E_Int));
  for (E_Int i = 0; i < M->nfaces; i++) {
    E_Int np = -1;
    E_Int *pn = cn->getFace(i, np, ngon, indPG);
    E_Int *ptr = &M->ngon[M->indPG[i]];
    for (E_Int j = 0; j < np; j++)
      ptr[j] = pn[j]-1;
  }

  M->x = (E_Float *)XMALLOC(M->npoints * sizeof(E_Float));
  M->y = (E_Float *)XMALLOC(M->npoints * sizeof(E_Float));
  M->z = (E_Float *)XMALLOC(M->npoints * sizeof(E_Float));
  memcpy(M->x, px, M->npoints*sizeof(E_Float));
  memcpy(M->y, py, M->npoints*sizeof(E_Float));
  memcpy(M->z, pz, M->npoints*sizeof(E_Float));

  // Parse comm data
  assert(PyList_Check(COMM));
  M->npatches = PyList_Size(COMM);
  M->patches = (Patch *)XMALLOC(M->npatches * sizeof(Patch));
  for (E_Int i = 0; i < M->npatches; i++) {
    PyObject *PATCH = PyList_GetItem(COMM, i);

    M->patches[i].nei_proc = PyLong_AsLong(PyList_GetItem(PATCH, 0));
    K_NUMPY::getFromNumpyArray(PyList_GetItem(PATCH, 1),
      M->patches[i].faces, M->patches[i].nfaces, false);
    
    // Zero-based
    for (E_Int j = 0; j < M->patches[i].nfaces; j++)
      M->patches[i].faces[j] -= 1;
    
    // Send buffer
    M->patches[i].sbuf_d = (E_Float *)XCALLOC(M->patches[i].nfaces,
      sizeof(E_Float));
  }

  // Adaptation trees
  M->cellTree = (Element **)XMALLOC(M->ncells * sizeof(Element *));
  M->faceTree = (Element **)XMALLOC(M->nfaces * sizeof(Element *));
  for (E_Int i = 0; i < M->ncells; i++) {
    M->cellTree[i] = (Element *)XMALLOC(sizeof(Element));
    Element *Elem = M->cellTree[i];
    Elem->children = NULL;
    Elem->nchildren = 0;
    Elem->parent = i;
    Elem->position = 0;
    Elem->type = -1;
    Elem->level = 0;
  }
  for (E_Int i = 0; i < M->nfaces; i++) {
    M->faceTree[i] = (Element *)XMALLOC(sizeof(Element));
    Element *Elem = M->faceTree[i];
    Elem->children = NULL;
    Elem->nchildren = 0;
    Elem->parent = i;
    Elem->position = 0;
    Elem->type = -1;
    Elem->level = 0;
  }

  M->ref_Tr = Tr;
  M->unref_Tr = Tr/2.5;

  // Set face types
  set_faces_type(M);

  // Set cell types
  set_cells_type(M);

  // Reorder cells
  reorder_cells(M);

  // Check everything is in place
  check_canon_cells(M);

  // Make hook
  // TODO(Imad): Python2
  PyObject *hook = PyCapsule_New((void *)M, "AMesh", NULL);

  RELEASESHAREDU(ARRAY, f, cn);

  return hook;
}

static
void get_face_leaves(Element **faceTree, E_Int face, std::vector<E_Int> &leaves)
{
  E_Int nchildren = faceTree[face]->nchildren;
  if (nchildren == 0) {
    leaves.push_back(face);
    return;
  }

  E_Int *children = faceTree[face]->children;
  for (E_Int i = 0; i < nchildren; i++)
    get_face_leaves(faceTree, children[i], leaves);
}

PyObject *K_XCORE::AdaptMesh(PyObject *self, PyObject *args)
{
  PyObject *AMESH, *FIELDS;
  
  if (!PYPARSETUPLE_(args, OO_, &AMESH, &FIELDS)) {
    RAISE("Wrong input.");
    return NULL;
  }

  // Unpack AMesh
  if (!PyCapsule_IsValid(AMESH, "AMesh")) {
    RAISE("Bad AMesh hook");
    return NULL;
  }

  AMesh *M = (AMesh *)PyCapsule_GetPointer(AMESH, "AMesh");

  // Parse fields
  assert(PyList_Check(FIELDS));
  E_Int nfields = PyList_Size(FIELDS);
  E_Float **fields = (E_Float **)XMALLOC(nfields * sizeof(E_Float *));
  for (E_Int i = 0; i < nfields; i++) {
    PyObject *FIELD = PyList_GetItem(FIELDS, i);
    E_Int nfld, nc;
    E_Int ret = K_NUMPY::getFromNumpyArray(FIELD, fields[i], nc, nfld, true);
    assert(ret == 1 && nfld == 1 && nc == M->ncells);
  }

  make_cell_centers(M);
  exchange_proc_data_d(M, M->cx, &M->xn);
  exchange_proc_data_d(M, M->cy, &M->yn);
  exchange_proc_data_d(M, M->cz, &M->zn);
  
  // Compute cell refinement levels
  compute_ref_data(M, fields, nfields);
  
  XFREE(fields);

  // Redistribute
  //AMesh *rM = redistribute_mesh(M);

  // Isolate refinement cells and faces
  std::vector<E_Int> ref_faces, ref_cells;

  // Resize structures for refinement
  get_ref_cells_and_faces(M, ref_cells, ref_faces);
  E_Int nref_cells = ref_cells.size();
  E_Int nref_faces = ref_faces.size();
  resize_data_for_refinement(M, nref_cells, nref_faces);
  printf("Refined cells: %d\n", nref_cells);

  /*
  ref_faces.clear();
  ref_cells.clear();
  for (E_Int i = 0; i < M->nfaces; i++) ref_faces.push_back(i);
  ref_cells.push_back(0);
  E_Int nref_cells = ref_cells.size();
  E_Int nref_faces = ref_faces.size();
  resize_data_for_refinement(M, nref_cells, nref_faces);
  */

  refine_faces(ref_faces, M);
  refine_cells(ref_cells, M);

  printf("Leaves: %d\n", M->ncells);

  std::map<E_Int, E_Int> epoints, efaces;

  std::vector<E_Int> ecells;
  for (E_Int i = 0; i < M->ncells; i++) {
    if (M->cellTree[i]->nchildren == 0)
      ecells.push_back(i);
  }

  E_Int necells = ecells.size();

  std::vector<E_Int> XCELLS(1, 0);
  std::vector<E_Int> XFACES(1, 0);

  std::vector<E_Int> NFACE, NGON;

  E_Int nefaces = 0;
  E_Int nepoints = 0;

  for (E_Int i = 0; i < necells; i++) {
    E_Int cell = ecells[i];
    E_Int stride = 0;
    for (E_Int j = M->indPH[cell]; j < M->indPH[cell+1]; j++) {
      E_Int face = M->nface[j];
      std::vector<E_Int> leaves;
      get_face_leaves(M->faceTree, face, leaves);
      assert(leaves.size() == 1 || leaves.size() == 4);
      stride += leaves.size();
      for (auto leaf : leaves)
        NFACE.push_back(leaf);
    }
    XCELLS.push_back(stride);
  }

  for (E_Int i = 0; i < necells; i++) XCELLS[i+1] += XCELLS[i];

  for (E_Int i = 0; i < necells; i++) {
    for (E_Int j = XCELLS[i]; j < XCELLS[i+1]; j++) {
      E_Int face = NFACE[j];
      if (efaces.find(face) == efaces.end()) {
        efaces[face] = nefaces++;
        for (E_Int k = M->indPG[face]; k < M->indPG[face+1]; k++) {
          E_Int point = M->ngon[k];
          NGON.push_back(point);
        }
        XFACES.push_back(M->indPG[face+1] - M->indPG[face]);
      }
    }
  }

  for (E_Int i = 0; i < nefaces; i++) XFACES[i+1] += XFACES[i];
  
  for (E_Int i = 0; i < nefaces; i++) {
    for (E_Int j = XFACES[i]; j < XFACES[i+1]; j++) {
      E_Int point = NGON[j];
      if (epoints.find(point) == epoints.end())
        epoints[point] = nepoints++;
    }
  }

  // Export
  const char *varStringOut = "CoordinateX,CoordinateY,CoordinateZ";

  PyObject *m = K_ARRAY::buildArray3(3, varStringOut, nepoints, necells,
    nefaces, "NGON", XFACES[nefaces], XCELLS[necells], 3, false, 3);

  FldArrayF *fo;
  FldArrayI *cno;
  K_ARRAY::getFromArray3(m, fo, cno);

  E_Int *ngono  = cno->getNGon();
  E_Int *nfaceo = cno->getNFace();
  E_Int *indPGo = cno->getIndPG();
  E_Int *indPHo = cno->getIndPH();

  E_Float *px = fo->begin(1);
  E_Float *py = fo->begin(2);
  E_Float *pz = fo->begin(3);
  for (auto pt: epoints) {
    E_Int lp = pt.second;
    assert(lp < nepoints);
    px[lp] = M->x[pt.first];
    py[lp] = M->y[pt.first];
    pz[lp] = M->z[pt.first];
  }

  for (E_Int i = 0; i < nefaces+1; i++) indPGo[i] = XFACES[i];
  for (E_Int i = 0; i < necells+1; i++) indPHo[i] = XCELLS[i];

  E_Int *ptr = ngono;

  for (E_Int i = 0; i < nefaces; i++) {
    E_Int start = XFACES[i];
    E_Int end = XFACES[i+1];
    for (E_Int j = start; j < end; j++) {
      *ptr = epoints[NGON[j]]+1;
      ptr++;
    }
  }
  
  ptr = nfaceo;

  for (E_Int i = 0; i < necells; i++) {
    E_Int start = XCELLS[i];
    E_Int end = XCELLS[i+1];
    for (E_Int j = start; j < end; j++) {
      *ptr = efaces[NFACE[j]]+1;
      ptr++;
    }
  }

  return m;
}
