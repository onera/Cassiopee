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

PyObject *K_XCORE::CreateAdaptMesh(PyObject *self, PyObject *args)
{
  PyObject *ARRAY, *OWN, *NEI, *COMM, *BCS, *NORMAL_2D;
  PyObject *GCELLS, *GFACES, *GPOINTS;
  E_Float Tr, Tu, eps, hmin, hmax;
  E_Int unref;

  if (!PYPARSETUPLE_(args, OOOO_ O_ RRRR_ R_ I_ O_ OOO_,
    &ARRAY, &OWN, &NEI, &COMM, &BCS,
    &Tr, &Tu, &eps, &hmin, &hmax, &unref, &NORMAL_2D,
    &GCELLS, &GFACES, &GPOINTS)) {
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

  // Parse parent elements
  E_Int size, nfld;
  ret = K_NUMPY::getFromNumpyArray(OWN, M->owner, size, nfld, false);
  assert(ret == 1 && size == M->nfaces && nfld == 1);
  ret = K_NUMPY::getFromNumpyArray(NEI, M->neigh, size, nfld, false);
  assert(ret == 1 && size == M->nfaces && nfld == 1);
 
  // Parse boundary conditions
  M->nbc = PyList_Size(BCS);
  //printf("%d boundary conditions found:\n", M->nbc);
  M->ptlists = (E_Int **)XCALLOC(M->nbc, sizeof(E_Int *));
  M->bcsizes = (E_Int *)XCALLOC(M->nbc, sizeof(E_Int));
  M->bcnames = (char **)XMALLOC(M->nbc * sizeof(char *));
  M->nbf = 0;
  for (E_Int i = 0; i < M->nbc; i++) {
    PyObject *BC = PyList_GetItem(BCS, i);
    PyObject *PTLIST = PyList_GetItem(BC, 0);
    PyObject *BCNAME = PyList_GetItem(BC, 1);
    K_NUMPY::getFromNumpyArray(PTLIST, M->ptlists[i], M->bcsizes[i], false);

    M->nbf += M->bcsizes[i];

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
    //printf("    %s\n", M->bcnames[i]);
  }

  // Parse comm data
  assert(PyList_Check(COMM));
  M->npatches = PyList_Size(COMM);
  M->patches = (Patch *)XMALLOC(M->npatches * sizeof(Patch));
  M->npf = 0;
  for (E_Int i = 0; i < M->npatches; i++) {
    PyObject *PATCH = PyList_GetItem(COMM, i);

    assert(PyList_Size(PATCH) == 3);

    M->patches[i].nei = PyLong_AsLong(PyList_GetItem(PATCH, 0));
    
    K_NUMPY::getFromNumpyArray(PyList_GetItem(PATCH, 1),
      M->patches[i].pf, M->patches[i].nf, false);
    
    K_NUMPY::getFromNumpyArray(PyList_GetItem(PATCH, 2),
      M->patches[i].pn, M->patches[i].nf, false);
    
    // Zero-based
    for (E_Int j = 0; j < M->patches[i].nf; j++)
      M->patches[i].pf[j] -= 1;
    
    M->npf += M->patches[i].nf;
    
    // Send buffers allocated on-demand
    M->patches[i].sbuf_i = NULL;
    M->patches[i].rbuf_i = NULL;
    M->patches[i].sbuf_f = NULL;
    M->patches[i].rbuf_f = NULL;
  }

  M->nif = M->nfaces - M->nbf - M->npf;

  // Interface buffers
  M->px   = (E_Float **)XCALLOC(M->npatches, sizeof(E_Float *));
  M->py   = (E_Float **)XCALLOC(M->npatches, sizeof(E_Float *));
  M->pz   = (E_Float **)XCALLOC(M->npatches, sizeof(E_Float *));
  M->pfld = (E_Float **)XCALLOC(M->npatches, sizeof(E_Float *));
  M->pref = (E_Int **)  XCALLOC(M->npatches, sizeof(E_Int *));
  M->plvl = (E_Int **)  XCALLOC(M->npatches, sizeof(E_Int *));

  // Adaptation trees
  M->cellTree = new Tree(M->ncells);
  M->faceTree = new Tree(M->nfaces);

  //init_mesh_numbering(M);

  M->Tr = Tr;
  M->Tu = Tu;
  M->eps = eps;
  M->hmin = hmin;
  M->hmax = hmax;
  M->unrefine = unref;
  
  // Normal to 2D plane
  if (NORMAL_2D != Py_None) {
    K_NUMPY::getFromNumpyArray(NORMAL_2D, M->mode_2D, size, nfld, false);
    assert(size == 3 && nfld == 1);
    //if (M->mode_2D) {
    //  printf("2D mode normal to (%f %f %f)\n", M->mode_2D[0], M->mode_2D[1],
    //    M->mode_2D[2]);
    //}
  }

  // Parse global cells / faces / points
  if (GCELLS != Py_None) {
    ret = K_NUMPY::getFromNumpyArray(GCELLS, M->gcells, size, nfld, false);
    assert(ret == 1 && size == M->ncells && nfld == 1);
    for (E_Int i = 0; i < M->ncells; i++) {
      M->CT->insert({M->gcells[i], i});
    }
  }

  if (GFACES != Py_None) {
    ret = K_NUMPY::getFromNumpyArray(GFACES, M->gfaces, size, nfld, false);
    assert(ret == 1 && size == M->nfaces && nfld == 1);
    for (E_Int i = 0; i < M->nfaces; i++) {
      M->gfaces[i] -= 1;
      assert(M->gfaces[i] >= 0);
      M->FT->insert({M->gfaces[i], i});
    }
  }

  if (GPOINTS != Py_None) {
    ret = K_NUMPY::getFromNumpyArray(GPOINTS, M->gpoints, size, nfld, false);
    assert(ret == 1 && size == M->npoints && nfld == 1);
    for (E_Int i = 0; i < M->npoints; i++) {
      M->gpoints[i] -= 1;
      assert(M->gpoints[i] >= 0);
      M->PT->insert({M->gpoints[i], i});
    }
  }

  // Set face types
  set_faces_type(M);

  // Set cell types
  set_cells_type(M);

  // Reorder cells
  reorder_cells(M);

  // Check everything is in place
  assert(check_canon_cells(M));

  // Make hook
  // TODO(Imad): Python2
  PyObject *hook = PyCapsule_New((void *)M, "AMesh", NULL);

  RELEASESHAREDU(ARRAY, f, cn);

  return hook;
}

static
PyObject *export_mesh_and_bcs(AMesh *M, E_Int closed_mesh);

static
void get_face_and_cell_centers(AMesh *M, PyObject *FC, PyObject *CX,
  PyObject *CY, PyObject *CZ)
{
  E_Int fsize, csize, nfld;
  E_Int *ptr;

  M->fc = (E_Float *)XRESIZE(M->fc, 3*M->nfaces * sizeof(E_Float));
  M->cx = (E_Float *)XRESIZE(M->cx, M->ncells * sizeof(E_Float));
  M->cy = (E_Float *)XRESIZE(M->cy, M->ncells * sizeof(E_Float));
  M->cz = (E_Float *)XRESIZE(M->cz, M->ncells * sizeof(E_Float));

  K_NUMPY::getFromNumpyArray(FC, ptr, fsize, nfld, true);
  assert(fsize == 3*M->nfaces);
  memcpy(M->fc, ptr, 3*M->nfaces * sizeof(E_Float));
  
  K_NUMPY::getFromNumpyArray(CX, ptr, csize, nfld, true);
  assert(csize == M->ncells);
  memcpy(M->cx, ptr, M->ncells * sizeof(E_Float));
  
  K_NUMPY::getFromNumpyArray(CY, ptr, csize, nfld, true);
  assert(csize == M->ncells);
  memcpy(M->cy, ptr, M->ncells * sizeof(E_Float));
  
  K_NUMPY::getFromNumpyArray(CZ, ptr, csize, nfld, true);
  assert(csize == M->ncells);
  memcpy(M->cz, ptr, M->ncells * sizeof(E_Float));
}

static
int FEQ(E_Float a, E_Float b)
{
  return fabs(a-b) < 1e-6;
}

int same_xyz(E_Float *A, E_Float *B, E_Int n)
{
  int same = 1;
  for (E_Int i = 0; i < n && same; i++) {
    if (!FEQ(A[i], B[i])){
      same = 0;
    }
  }
  return same;
}

PyObject *K_XCORE::AdaptMesh(PyObject *self, PyObject *args)
{
  PyObject *AMESH, *FC, *CX, *CY, *CZ;
  
  if (!PYPARSETUPLE_(args, OOO_ OO_, &AMESH, &FC, &CX, &CY, &CZ)) {
    RAISE("Wrong input.");
    return NULL;
  }

  // Unpack AMesh
  if (!PyCapsule_IsValid(AMESH, "AMesh")) {
    RAISE("Bad AMesh hook");
    return NULL;
  }

  AMesh *M = (AMesh *)PyCapsule_GetPointer(AMESH, "AMesh");

  get_face_and_cell_centers(M, FC, CX, CY, CZ);


  //M = Redistribute_mesh(M);

  std::vector<E_Int> ref_faces, ref_cells;

  get_ref_faces_and_cells(M, ref_faces, ref_cells);
  
  M->prev_ncells = M->ncells;
  M->prev_nfaces = M->nfaces;
  M->prev_npoints = M->npoints;

  refine_mesh(M, ref_faces, ref_cells);

  update_patch_faces_after_ref(M);

  MPI_Barrier(MPI_COMM_WORLD);

  E_Int gncells;
  MPI_Allreduce(&M->ncells, &gncells, 1, XMPI_INT, MPI_SUM, MPI_COMM_WORLD);

  if (M->pid == 0)
    printf("total leaves: %d\n", gncells);

  for (E_Int i = 0; i < M->ncells; i++) M->cellTree->state_[i] = UNTOUCHED;
  for (E_Int i = 0; i < M->nfaces; i++) M->faceTree->state_[i] = UNTOUCHED;

  E_Int closed_mesh = 1;

  PyObject *out = export_mesh_and_bcs(M, closed_mesh);

  return out;
}

//#define PRINTPATCH

static
PyObject *export_mesh_and_bcs(AMesh *M, E_Int closed_mesh)
{
  // Count nface size
  E_Int nface_size = 0;

  M->closed_indPH = (E_Int *)XRESIZE(M->closed_indPH, (M->ncells+1)*sizeof(E_Int));
  M->closed_indPH[0] = 0;

  for (E_Int i = 0; i < M->ncells; i++) {
    E_Int nf = -1;
    E_Int pf[24];
    // TODO(Imad): a function that returns nf only
    get_full_cell(i, M, nf, pf);

    nface_size += nf;
    M->closed_indPH[i+1] = nf;
  }

  for (E_Int i = 0; i < M->ncells; i++)
    M->closed_indPH[i+1] += M->closed_indPH[i];

  assert(M->closed_indPH[M->ncells] == nface_size);

  E_Int ngon_size = -1;

  if (closed_mesh) {
    puts("Exporting closed mesh");
    close_mesh(M);
    ngon_size = M->closed_indPG[M->nfaces];
  } else {
    puts("Exporting non-closed mesh");
    ngon_size = M->indPG[M->nfaces];
  }

  //M->closed_nface = (E_Int *)XRESIZE(M->closed_nface, nface_size * sizeof(E_Int));

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
  
  /*
  indPHo[0] = 0;
  for (E_Int i = 0; i < M->ncells; i++) {
    E_Int nf = -1;
    E_Int pf[24];
    get_full_cell(i, M, nf, pf);
    indPHo[i+1] = nf;
  }
  for (E_Int i = 0; i < M->ncells; i++) indPHo[i+1] += indPHo[i];
  */
  memcpy(indPHo, M->closed_indPH, (M->ncells+1)*sizeof(E_Int));
  assert(indPHo[M->ncells] == nface_size);

  for (E_Int i = 0; i < M->ncells; i++) {
    E_Int nf = -1;
    E_Int pf[24];
    get_full_cell(i, M, nf, pf);
    E_Int *ptr = &nfaceo[indPHo[i]];
    //E_Int *p = &M->closed_nface[M->closed_indPH[i]];
    for (E_Int j = 0; j < nf; j++) {
      ptr[j] = pf[j]+1;
      //p[j] = pf[j];
    }
  }

  if (closed_mesh) {
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

  assert(K_CONNECT::check_open_cells(*cno, NULL) == 0);
  assert(K_CONNECT::check_overlapping_cells(*cno) == 0);

  // BCs
  PyObject *BCS = PyList_New(0);

  npy_intp dims[2];
  
  for (E_Int i = 0; i < M->nbc; i++) {
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

#ifdef PRINTPATCH
  // Patch faces as BCs
  for (E_Int i = 0; i < M->npatches; i++) {
    Patch *P = &M->patches[i];

    dims[0] = (npy_intp)P->nf;
    dims[1] = 1;

    PyArrayObject *PL = (PyArrayObject *)PyArray_SimpleNew(1, dims, E_NPY_INT);
    E_Int *op = P->pf;
    E_Int *np = (E_Int *)PyArray_DATA(PL);
    for (E_Int j = 0; j < P->nf; j++)
      *np++ = op[j] + 1;

    char name[64];
    sprintf(name, "proc%d-%d", M->pid, P->nei);

    PyObject *tpl = Py_BuildValue("[Os]", (PyObject *)PL, name);
    Py_DECREF(PL);
    PyList_Append(BCS, tpl);
    Py_DECREF(tpl);
  }
#endif

  // parent elements
  dims[0] = (npy_intp)M->nfaces;
  dims[1] = 1;

  PyArrayObject *OWN = (PyArrayObject *)PyArray_SimpleNew(1, dims, E_NPY_INT);
  PyArrayObject *NEI = (PyArrayObject *)PyArray_SimpleNew(1, dims, E_NPY_INT);
  E_Int *owner = (E_Int *)PyArray_DATA(OWN);
  E_Int *neigh = (E_Int *)PyArray_DATA(NEI);

  memcpy(owner, M->owner, M->nfaces*sizeof(E_Int));
  memcpy(neigh, M->neigh, M->nfaces*sizeof(E_Int));

  PyObject *out = PyList_New(0);
  PyList_Append(out, m);
  PyList_Append(out, BCS);
  PyList_Append(out, (PyObject *)OWN);
  PyList_Append(out, (PyObject *)NEI);
  Py_DECREF(m);
  Py_DECREF(BCS);
  Py_DECREF(OWN);
  Py_DECREF(NEI);

  return out;
}
