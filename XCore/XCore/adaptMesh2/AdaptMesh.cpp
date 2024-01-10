#include "Proto.h"

// 1. Make a mesh with bcs
// 2. Create AdaptMesh
// 3. Renumber mesh: internal faces - boundary faces

// 3. Convert to foam mesh
// 4. Map prev solution onto current mesh
// 5. Run for niters

// 6. Convert to cgns
// 7. Adapt (refine - unrefine)
// 8. Shrink adapt data
// 9. Renumber mesh: internal faces - boundary faces

// 10. Repeat steps 3 to 9 until end time

#define MINREF -10
#define RAISE(error) PyErr_SetString(PyExc_ValueError, (error))

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

PyObject *K_XCORE::CreateAdaptMesh(PyObject *self, PyObject *args)
{
  PyObject *ARRAY, *COMM, *BCS;
  E_Float Tr;
  if (!PYPARSETUPLE_(args, OOO_ R_, &ARRAY, &COMM, &BCS, &Tr)) {
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

  E_Float *px = f->begin(posx);
  E_Float *py = f->begin(posy);
  E_Float *pz = f->begin(posz);

  // Init mesh
  AMesh *M = init_mesh(*cn, px, py, pz, f->getSize());
 
  // Parse boundary conditions
  M->nbc = PyList_Size(BCS);
  printf("%d boundary conditions found:\n", M->nbc);
  M->ptlists = (E_Int **)XCALLOC(M->nbc, sizeof(E_Int *));
  M->bcsizes = (E_Int *)XCALLOC(M->nbc, sizeof(E_Int));
  M->bcnames = (char **)XMALLOC(M->nbc * sizeof(char *));
  for (E_Int i = 0; i < M->nbc; i++) {
    PyObject *BC = PyList_GetItem(BCS, i);
    PyObject *PTLIST = PyList_GetItem(BC, 0);
    PyObject *BCNAME = PyList_GetItem(BC, 1);
    K_NUMPY::getFromNumpyArray(PTLIST, M->ptlists[i], M->bcsizes[i], false);

    // Zero-based
    E_Int *ptr = M->ptlists[i];
    for (E_Int j = 0; j < M->bcsizes[i]; j++) ptr[j] -= 1;

    M->bcnames[i] = (char *)XMALLOC(128*sizeof(char));

#if PY_VERSION_HEX >= 0x03000000
    char *tmp = (char *)PyUnicode_AsUTF8(BCNAME);
#else
    char *tmp = PyString_AsString(BCNAME);
#endif

    strcpy(M->bcnames[i], tmp);
    printf("    %s\n", M->bcnames[i]);
  }

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

  // Init owner and neigh
  M->owner = (E_Int *)XMALLOC(M->nfaces * sizeof(E_Int));
  M->neigh = (E_Int *)XMALLOC(M->nfaces * sizeof(E_Int));
  memset(M->owner, -1, M->nfaces * sizeof(E_Int));
  memset(M->neigh, -1, M->nfaces * sizeof(E_Int));
  M->nif = 0;
  for (E_Int i = 0; i < M->ncells; i++) {
    for (E_Int j = M->indPH[i]; j < M->indPH[i+1]; j++) {
      E_Int face = M->nface[j];
      if (M->owner[face] == -1) {
        M->owner[face] = i;
      } else {
        M->neigh[face] = i;
        M->nif++;
      }
    }
  }

  // Adaptation trees
  M->cellTree = (Element **)XMALLOC(M->ncells * sizeof(Element *));
  M->faceTree = (Element **)XMALLOC(M->nfaces * sizeof(Element *));
  tree_init(M->cellTree, M->ncells);
  tree_init(M->faceTree, M->nfaces);

  // TODO(Imad): Check mesh validity
  /*
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
  */

  init_mesh_numbering(M);

  // Build owner and neigh
  ret = Orient_boundary(M);

  if (ret != 0) {
    RELEASESHAREDU(ARRAY, f, cn);
    // TODO(Imad): XFREE(M)
    RAISE("Failed boundary orientation.");
    return NULL;
  }

  Build_own_nei(M);

  M->ref_Tr = Tr;
  M->unref_Tr = Tr/2.5;

  // Renumber
  // TODO(Imad): optional
  //renumber_mesh(M);

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
const char *TYPE(E_Int type)
{
  switch (type) {
    case HEXA: return "HEXA";
    case TETRA: return "TETRA";
    case PENTA: return "PENTA";
    case PYRA: return "PYRA";
    default: assert(0);
  }
}

static
void print_elem(E_Int id, Element **tree)
{
  Element *E = tree[id];
  printf("%d [%s, %d]", id, TYPE(E->type), E->parent);
}

static
void print_elem_type(E_Int id, Element **tree)
{
  Element *E = tree[id];
  printf("[%s]", TYPE(E->type));
}

static
PyObject *export_mesh_and_bcs(AMesh *M, E_Int closed_mesh);

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
  //for (E_Int i = M->indPH[0]; i < M->indPH[1]; i++)
  //  ref_faces.push_back(M->nface[i]);
  //ref_cells.push_back(0);
  
  E_Int nref_cells = ref_cells.size();
  E_Int nref_faces = ref_faces.size();
  resize_data_for_refinement(M, nref_cells, nref_faces);
  printf("Refined cells: %d\n", nref_cells);

  refine_faces(ref_faces, M);
  refine_cells(ref_cells, M);

  printf("Leaves: %d\n", M->ncells);

  update_boundary_faces(M);

  init_mesh_numbering(M);

  //renumber_mesh(M);

  E_Int closed_mesh = 1;

  PyObject *out = export_mesh_and_bcs(M, closed_mesh);

  return out;
}

static
PyObject *export_mesh_and_bcs(AMesh *M, E_Int closed_mesh)
{
  // Count nface size
  E_Int nface_size = 0;

  for (E_Int i = 0; i < M->ncells; i++) {
    E_Int nf = -1;
    E_Int pf[24];
    // TODO(Imad): a function that returns nf only
    get_full_cell(i, M, nf, pf);
    nface_size += nf;
  }

  E_Int ngon_size = -1;

  if (closed_mesh) {
    puts("Exporting closed mesh");
    close_mesh(M);
    ngon_size = M->closed_indPG[M->nfaces];
  } else {
    puts("Exporting non-closed mesh");
    ngon_size = M->indPG[M->nfaces];
  }

  const char *varStringOut = "CoordinateX,CoordinateY,CoordinateZ";

  PyObject *m = K_ARRAY::buildArray3(3, varStringOut, M->npoints, M->ncells,
    M->nfaces, "NGON", ngon_size, nface_size, 3, false, 3);

  FldArrayF *fo;
  FldArrayI *cno;
  K_ARRAY::getFromArray3(m, fo, cno);

  E_Int *ngono  = cno->getNGon();
  E_Int *nfaceo = cno->getNFace();
  E_Int *indPGo = cno->getIndPG();
  E_Int *indPHo = cno->getIndPH();

  E_Int *ptr = nfaceo;

  for (E_Int i = 0; i < M->ncells; i++) {
    get_full_cell(i, M, indPHo[i+1], ptr);
    ptr += indPHo[i+1];
  }

  for (E_Int i = 0; i < nface_size; i++) nfaceo[i] += 1;

  indPHo[0] = 0;
  for (E_Int i = 0; i < M->ncells; i++) indPHo[i+1] += indPHo[i];

  assert(indPHo[M->ncells] == nface_size);

  if (close_mesh) {
    for (E_Int i = 0; i < ngon_size; i++)
      ngono[i] = M->closed_ngon[i] + 1;
    
    memcpy(indPGo, M->closed_indPG, (M->nfaces+1) * sizeof(E_Int));
  } else {
    for (E_Int i = 0; i < ngon_size; i++)
      ngono[i] = M->ngon[i] + 1;
    
    memcpy(indPGo, M->indPG, (M->nfaces+1) * sizeof(E_Int));
  }

  E_Float *px = fo->begin(1);
  E_Float *py = fo->begin(2);
  E_Float *pz = fo->begin(3);
  memcpy(px, M->x, M->npoints * sizeof(E_Float));
  memcpy(py, M->y, M->npoints * sizeof(E_Float));
  memcpy(pz, M->z, M->npoints * sizeof(E_Float));

  // BCs
  PyObject *BCS = PyList_New(0);

  for (E_Int i = 0; i < M->nbc; i++) {
    npy_intp dims[2];
    dims[0] = (npy_intp)M->bcsizes[i];
    dims[1] = 1;

    PyArrayObject *PL = (PyArrayObject *)PyArray_SimpleNew(1, dims, E_NPY_INT);
    E_Int *op = M->ptlists[i];
    E_Int *np = (E_Int *)PyArray_DATA(PL);
    for (E_Int j = 0; j < M->bcsizes[i]; j++)
      *np++ = op[j] + 1;

    PyObject *tpl = Py_BuildValue("[Os]", (PyObject *)PL, M->bcnames[i]);
    Py_DECREF(PL);
    PyList_Append(BCS, tpl);
    Py_DECREF(tpl);
  }

  PyObject *out = PyList_New(0);
  PyList_Append(out, m);
  PyList_Append(out, BCS);
  Py_DECREF(m);
  Py_DECREF(BCS);

  return out;
}