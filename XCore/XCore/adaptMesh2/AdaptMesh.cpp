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
  PyObject *ARRAY, *OWN, *NEI, *COMM, *BCS;
  E_Float Tr, Tu;
  if (!PYPARSETUPLE_(args, OOO_ OO_ RR_, &ARRAY, &OWN, &NEI, &COMM, &BCS, &Tr, &Tu)) {
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
  printf("%d boundary conditions found:\n", M->nbc);
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
    printf("    %s\n", M->bcnames[i]);
  }
  M->nif = M->nfaces - M->nbf;

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
  M->cellTree = new Tree(M->ncells);
  M->faceTree = new Tree(M->nfaces);

  //printf("Boundary faces: %d\n", M->nbf);
  //printf("Internal faces: %d\n", M->nif);
  //init_mesh_numbering(M);

  M->Tr = Tr;
  M->Tu = Tu;

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

  printf("    Initial cells: %d\n", M->ncells);

  // Isolate refinement cells and faces
  std::vector<E_Int> ref_faces, ref_cells;
  std::vector<E_Int> unref_faces, unref_cells;

  get_ref_cells_and_faces(M, ref_cells, ref_faces, unref_cells, unref_faces);

  //for (E_Int i = 0; i < M->nfaces; i++) M->faceTree->state_[i] = UNTOUCHED;
  //for (E_Int i = 0; i < M->ncells; i++) M->cellTree->state_[i] = UNTOUCHED;

  //for (E_Int i = M->indPH[0]; i < M->indPH[1]; i++)
  //  ref_faces.push_back(M->nface[i]);
  //ref_cells.push_back(0);

  size_t nref_cells = ref_cells.size();
  size_t nref_faces = ref_faces.size();
  size_t nunref_cells = unref_cells.size();
  size_t nunref_faces = unref_faces.size();

  printf("    Cells tagged for refinement: %d\n", (E_Int)nref_cells);
  printf("    Cells tagged for unrefinement: %d\n", (E_Int)nunref_cells);
  
  if (nunref_cells > 0) {

    // Sort unref cells/faces by decreasing level
    std::sort(unref_faces.begin(), unref_faces.end(), [&] (E_Int i, E_Int j)
    {
      return M->faceTree->level(i) > M->faceTree->level(j);
    });

    std::sort(unref_cells.begin(), unref_cells.end(), [&] (E_Int i, E_Int j)
    {
      return M->cellTree->level(i) > M->cellTree->level(j);
    });

    std::vector<E_Int> face_distrib(1, 0), cell_distrib(1, 0);
    std::vector<E_Int> face_lvl, cell_lvl;

    E_Int max_lvl = M->faceTree->level(unref_faces[0]);
    face_lvl.push_back(max_lvl);
    for (size_t i = 1; i < unref_faces.size(); i++) {
      if (M->faceTree->level(unref_faces[i]) < max_lvl) {
        max_lvl = M->faceTree->level(unref_faces[i]);
        face_distrib.push_back(i);
        face_lvl.push_back(max_lvl);
      }
    }
    face_distrib.push_back(unref_faces.size());

    max_lvl = M->cellTree->level(unref_cells[0]);
    cell_lvl.push_back(max_lvl);
    for (size_t i = 1; i < unref_cells.size(); i++) {
      if (M->cellTree->level(unref_cells[i]) < max_lvl) {
        max_lvl = M->cellTree->level(unref_cells[i]);
        cell_distrib.push_back(i);
        cell_lvl.push_back(max_lvl);
      }
    }
    cell_distrib.push_back(unref_cells.size());

    assert(face_distrib.size() == cell_distrib.size());
    
    /*
    size_t cell_start = 0;
    while (cell_lvl[cell_start] > face_lvl[0]) {
      unrefine_cells(unref_cells, cell_distrib[cell_start],
        cell_distrib[cell_start+1], M);
      cell_start++;
    }
    */

    for (size_t i = 0; i < face_distrib.size()-1; i++) {
      unrefine_cells(unref_cells, cell_distrib[i], cell_distrib[i+1], M);
      unrefine_faces(unref_faces, face_distrib[i], face_distrib[i+1], M);
      
      for (E_Int j = cell_distrib[i]; j < cell_distrib[i+1]; j++) {
        E_Int ucell = unref_cells[j];

        Children *cur = M->cellTree->children(ucell);
        Children *prev = cur->next;
        M->cellTree->children_[ucell] = prev;
  
        // free
        XFREE(cur);
      }
      
      for (E_Int j = face_distrib[i]; j < face_distrib[i+1]; j++) {
        E_Int uface = unref_faces[j];

        Children *cur = M->faceTree->children(uface);
        if (cur == NULL) {
          assert(M->owner[uface] == -1);
          assert(M->neigh[uface] == -1);
          continue;
        }
        Children *prev = cur->next;
        M->faceTree->children_[uface] = prev;
  
        // free
        XFREE(cur);
      }
    }

    update_boundary_faces(M);

    // Renumber cells
    E_Int nc = 0;
    E_Int sizeNFace = 0;
    Tree *CT = M->cellTree;
    std::vector<E_Int> new_cells(M->ncells, -1);
    for (E_Int i = 0; i < M->ncells; i++) {
      if (CT->state(i) != GONE) {
        new_cells[i] = nc++;
        sizeNFace += get_stride(i, M->indPH);
      }
    }

    // Renumber faces
    E_Int nf = 0;
    E_Int sizeNGon = 0;
    Tree *FT = M->faceTree;
    std::vector<E_Int> new_faces(M->nfaces, -1);
    for (E_Int i = 0; i < M->nfaces; i++) {
      if (FT->state(i) != GONE) {
        new_faces[i] = nf++;
        sizeNGon += get_stride(i, M->indPG);
      }
    }

    renumber_mesh(M, new_cells, new_faces, nc, nf, sizeNFace, sizeNGon);

    CT->compress(new_cells, nc);
    
    FT->compress(new_faces, nf);

    std::vector<E_Int> tmp_cells(ref_cells);

    // Renumber ref cells/faces
    for (size_t i = 0; i < ref_cells.size(); i++) {
      E_Int cell = tmp_cells[i];
      ref_cells[i] = new_cells[cell];
      assert(ref_cells[i] != -1);
      assert(M->cellTree->state(ref_cells[i]) == UNTOUCHED);
    }

    std::vector<E_Int> tmp_faces(ref_faces);

    for (size_t i = 0; i < ref_faces.size(); i++) {
      E_Int face = tmp_faces[i];
      ref_faces[i] = new_faces[face];
      assert(ref_faces[i] != -1);
      assert(M->faceTree->state(ref_faces[i]) == UNTOUCHED);
    }

    printf("    Cells after unrefinement: %d\n", M->ncells);

    //goto finish;
  }
 
 
  if (nref_cells > 0) {

    resize_data_for_refinement(M, nref_cells, nref_faces);

    // Sort ref cells/faces by increasing level
    std::sort(ref_faces.begin(), ref_faces.end(), [&] (E_Int i, E_Int j)
    {
      return M->faceTree->level(i) < M->faceTree->level(j);
    });

    std::sort(ref_cells.begin(), ref_cells.end(), [&] (E_Int i, E_Int j)
    {
      return M->cellTree->level(i) < M->cellTree->level(j);
    });

    std::vector<E_Int> face_distrib(1, 0), cell_distrib(1, 0);
    std::vector<E_Int> face_lvl, cell_lvl;

    E_Int min_lvl = M->faceTree->level(ref_faces[0]);
    face_lvl.push_back(min_lvl);
    for (size_t i = 1; i < ref_faces.size(); i++) {
      if (M->faceTree->level(ref_faces[i]) > min_lvl) {
        min_lvl = M->faceTree->level(ref_faces[i]);
        face_distrib.push_back(i);
        face_lvl.push_back(min_lvl);
      }
    }
    face_distrib.push_back(ref_faces.size());

    min_lvl = M->cellTree->level(ref_cells[0]);
    cell_lvl.push_back(min_lvl);
    for (size_t i = 1; i < ref_cells.size(); i++) {
      if (M->cellTree->level(ref_cells[i]) > min_lvl) {
        min_lvl = M->cellTree->level(ref_cells[i]);
        cell_distrib.push_back(i);
        cell_lvl.push_back(min_lvl);
      }
    }
    cell_distrib.push_back(ref_cells.size());

    size_t cell_start = 0;
    while (cell_lvl[cell_start] < face_lvl[0]) {
      refine_cells(ref_cells, cell_distrib[cell_start],
        cell_distrib[cell_start+1], M);
      cell_start++;
    }

    for (size_t i = 0; i < face_distrib.size()-1; i++) {
      refine_faces(ref_faces, face_distrib[i], face_distrib[i+1], M);
      refine_cells(ref_cells, cell_distrib[cell_start],
        cell_distrib[cell_start+1], M);
      cell_start++;
    }

    assert(cell_start == cell_distrib.size()-1);

    printf("    Cells after refinement: %d\n", M->ncells);

    update_boundary_faces(M);

  }

  init_mesh_numbering(M);

  //finish:

  for (E_Int i = 0; i < M->ncells; i++) M->cellTree->state_[i] = UNTOUCHED;
  for (E_Int i = 0; i < M->nfaces; i++) M->faceTree->state_[i] = UNTOUCHED;

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

  indPHo[0] = 0;
  for (E_Int i = 0; i < M->ncells; i++) {
    E_Int nf = -1;
    E_Int pf[24];
    get_full_cell(i, M, nf, pf);
    indPHo[i+1] = nf;
  }
  for (E_Int i = 0; i < M->ncells; i++) indPHo[i+1] += indPHo[i];
  assert(indPHo[M->ncells] == nface_size);

  for (E_Int i = 0; i < M->ncells; i++) {
    E_Int nf = -1;
    E_Int pf[24];
    get_full_cell(i, M, nf, pf);
    E_Int *ptr = &nfaceo[indPHo[i]];
    for (E_Int j = 0; j < nf; j++)
      ptr[j] = pf[j]+1;
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
