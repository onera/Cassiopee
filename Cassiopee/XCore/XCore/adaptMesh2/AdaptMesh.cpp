/*    
    Copyright 2013-2024 Onera.

    This file is part of Cassiopee.

    Cassiopee is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    Cassiopee is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Cassiopee.  If not, see <http://www.gnu.org/licenses/>.
*/
#include "Proto.h"

PyObject *K_XCORE::AdaptMesh(PyObject *self, PyObject *args)
{
  PyObject *AMESH;
  
  if (!PYPARSETUPLE_(args, O_, &AMESH)) {
    RAISE("Wrong input.");
    return NULL;
  }

  // Unpack AMesh
  if (!PyCapsule_IsValid(AMESH, "AMesh")) {
    RAISE("Bad AMesh hook");
    return NULL;
  }

  AMesh *M = (AMesh *)PyCapsule_GetPointer(AMESH, "AMesh");

  for (E_Int i = 0; i < M->ncells; i++) M->cellTree->state_[i] = UNTOUCHED;
  for (E_Int i = 0; i < M->nfaces; i++) M->faceTree->state_[i] = UNTOUCHED;

  //M = Redistribute_mesh(M);

  std::vector<E_Int> ref_faces, ref_cells;

  /*
  ref_cells.push_back(0);
  E_Int *pf = get_facets(0, M->nface, M->indPH);
  for (E_Int i = 0; i < 6; i++) ref_faces.push_back(pf[i]);
  */


  get_ref_faces_and_cells(M, ref_faces, ref_cells);

  printf("%d -> nref_cells: %zu - nref_faces: %zu\n", M->pid, ref_cells.size(),
    ref_faces.size());
  
  M->prev_ncells = M->ncells;
  M->prev_nfaces = M->nfaces;
  M->prev_npoints = M->npoints;

  resize_data_for_refinement(M, ref_cells.size(), ref_faces.size());

  if (ref_faces.size()) {
    std::sort(ref_faces.begin(), ref_faces.end(), [&] (E_Int i, E_Int j) {
      return M->faceTree->level(i) < M->faceTree->level(j);
    });
  }

  if (ref_cells.size()) {
    std::sort(ref_cells.begin(), ref_cells.end(), [&] (E_Int i, E_Int j) {
      return M->cellTree->level(i) < M->cellTree->level(j);
    });
  }

  compute_ref_cells_centers(M, ref_cells);

  refine_mesh(M, ref_faces, ref_cells);

  update_patch_faces_after_ref(M);

  E_Int gncells;
  MPI_Allreduce(&M->ncells, &gncells, 1, XMPI_INT, MPI_SUM, MPI_COMM_WORLD);

  if (M->pid == 0)
    printf("Total leaves: " SF_D_ "\n", gncells);

  return Py_None;
}

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

  M->onc = M->ncells;
  M->onf = M->nfaces;
  M->onp = M->npoints;

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
