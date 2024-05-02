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
#include "proto.h"

static
void shift_data(mesh *M)
{
  for (E_Int i = 0; i < M->nfaces; i++) {
    for (E_Int j = M->xfaces[i]; j < M->xfaces[i+1]; j++)
      M->NGON[j] -= 1;
  }

  for (E_Int i = 0; i < M->ncells; i++) {
    for (E_Int j = M->xcells[i]; j < M->xcells[i+1]; j++)
      M->NFACE[j] -= 1;
  }
}

#define MANDATORY 0
#define OPTIONAL 1

/*
  Returns 0 if:
   - mandatory key found in dict
   - optional key not found in dict
  Returns 1 if:
   - unknown key mode (not mandatory and not optional)
   - mandatory key not found in dict
*/
static
E_Int parse_dict(PyObject *dict, const char *key, E_Int mode,
  PyObject **obj)
{
  if (mode != MANDATORY && mode != OPTIONAL) {
    PyErr_Format(PyExc_ValueError,
      "adaptMesh(): unknown mode for key %s.", key);
    return 1;
  }

  PyObject *strObj;
#if PY_VERSION_HEX >= 0x03000000
  strObj = PyUnicode_FromString(key);
#else
  strObj = PyString_FromString(key);
#endif

  E_Int ret = 0;
  *obj = PyDict_GetItem(dict, strObj);
  if (*obj == NULL && mode == MANDATORY) {
    PyErr_Format(PyExc_KeyError, "Missing mandatory key %s in dictionary", key);
    ret = 1;
  }

  Py_DECREF(strObj);
  return ret;
}

PyObject *K_XCORE::adaptMesh(PyObject *self, PyObject *args)
{
  PyObject *arr, *comm_data, *solc, *gcells, *gfaces, *gpoints, *adaptDict;

  if (!PyArg_ParseTuple(args, "OOOOOOO", &arr,
    &comm_data, &solc, &gcells, &gfaces, &gpoints, &adaptDict)) {
    PyErr_SetString(PyExc_ValueError, "adaptMesh(): wrong input.");
    return NULL;
  }

  if (!PyDict_Check(adaptDict)) {
    PyErr_SetString(PyExc_ValueError,
      "adaptMesh(): adaptDict should be a python dictionary.");
      return NULL;
  }

  // Input parameters
  E_Int   Gmax;           // Mandatory
  E_Float Tr;             // Mandatory
  E_Int   conformize;     // Optional
  E_Int   sensor;         // Optional
  E_Int   iso_mode;       // Optional
  E_Float *freezeVector;  // Optional

  // Default values for optional parameters
  conformize   = 1;     // Conformal final mesh
  sensor       = 0;     // Metric sensor
  iso_mode     = 0;     // Directional adaptation
  freezeVector = NULL;  // Allow refinement in every direction

  // Parse dictionary
  PyObject *obj = NULL;
  E_Int ret;

  // Gmax
  ret = parse_dict(adaptDict, "Gmax", MANDATORY, &obj);
  if (ret) return NULL;
  Gmax = PyLong_AsLong(obj);

  // Tr
  ret = parse_dict(adaptDict, "Tr", MANDATORY, &obj);
  if (ret) return NULL;
  Tr = PyFloat_AsDouble(obj);

  // iso_mode
  parse_dict(adaptDict, "iso_mode", OPTIONAL, &obj);
  if (obj) iso_mode = PyLong_AsLong(obj);

  // conformize
  parse_dict(adaptDict, "conformize", OPTIONAL, &obj);
  if (obj) conformize = PyLong_AsLong(obj);

  // sensor
  parse_dict(adaptDict, "sensor", OPTIONAL, &obj);
  if (obj) sensor = PyLong_AsLong(obj);

  // freeze vector
  parse_dict(adaptDict, "freezeVector", OPTIONAL, &obj);
  if (obj) {
    E_Int size, nfld, ret;
    ret = K_NUMPY::getFromNumpyArray(obj, freezeVector, size, nfld, true);
    if (ret == 0 || nfld != 1 || size != 3) {
      PyErr_SetString(PyExc_ValueError, "adaptMesh(): bad freezeVector.");
      return NULL;
    }
  }

  // Check array
  E_Int ni, nj, nk;
  K_FLD::FldArrayF* f; K_FLD::FldArrayI* cn;
  char* varString; char* eltType;
  ret = K_ARRAY::getFromArray3(arr, varString, f, ni, nj, nk, cn, eltType);
  
  if (ret <= 0) {
    PyErr_SetString(PyExc_TypeError, "adaptMesh(): only for NGons.");
    return NULL;
  }

  if (ret == 1) { 
    PyErr_SetString(PyExc_TypeError, "adaptMesh(): only for NGons."); 
    RELEASESHAREDS(arr, f);
    return NULL; 
  }

  E_Int posx = K_ARRAY::isCoordinateXPresent(varString);
  E_Int posy = K_ARRAY::isCoordinateYPresent(varString);
  E_Int posz = K_ARRAY::isCoordinateZPresent(varString);
  if (posx == -1 || posy == -1 || posz == -1)
  {
    RELEASESHAREDB(ret, arr, f, cn);
    PyErr_SetString(PyExc_ValueError,
                    "adaptMesh(): can't find coordinates in array.");
    return NULL;
  }
  posx++; posy++; posz++;

  // Init mesh
  mesh *M = new mesh();
  M->NGON = cn->getNGon();
  M->NFACE = cn->getNFace();
  M->xfaces = cn->getIndPG();
  M->xcells = cn->getIndPH();
  M->Gmax = Gmax;
  M->Tr = Tr;
  M->iso_mode = iso_mode;
  M->sensor = sensor;

  // TODO(Imad): avoid copying
  // Coordinates
  M->npoints = f->getSize();
  M->xyz = (E_Float *)XCALLOC(3*M->npoints, sizeof(E_Float));
  E_Float *X = f->begin(1);
  E_Float *Y = f->begin(2);
  E_Float *Z = f->begin(3);
  for (E_Int i = 0; i < M->npoints; i++) {
    E_Float *px = &M->xyz[3*i];
    px[0] = X[i];
    px[1] = Y[i];
    px[2] = Z[i];
  }

  // Parse my global cells, faces and points
  E_Int nfld;
  K_NUMPY::getFromNumpyArray(gcells, M->gcells, M->ncells, nfld, true);
  for (E_Int i = 0; i < M->ncells; i++) {
    M->CT[M->gcells[i]] = i;
  }

  K_NUMPY::getFromNumpyArray(gfaces, M->gfaces, M->nfaces, nfld, true);
  for (E_Int i = 0; i < M->nfaces; i++) {
    M->FT[M->gfaces[i]] = i;
  }

  K_NUMPY::getFromNumpyArray(gpoints, M->gpoints, M->npoints, nfld, true);
  for (E_Int i = 0; i < M->npoints; i++) {
    M->PT[M->gpoints[i]] = i;
  }

  E_Int gncells = 0;
  MPI_Allreduce(&M->ncells, &gncells, 1, XMPI_INT, MPI_SUM, MPI_COMM_WORLD);

  if (M->pid == 0)
    printf("Total number of cells: " SF_D_ "\n", gncells);

  shift_data(M);

  // Init parent elements
  M->owner = (E_Int *)XCALLOC(M->nfaces, sizeof(E_Int));
  M->neigh = (E_Int *)XCALLOC(M->nfaces, sizeof(E_Int));
  for (E_Int i = 0; i < M->nfaces; i++) {
    M->owner[i] = -1;
    M->neigh[i] = -1;
  }
  for (E_Int i = 0; i < M->ncells; i++) {
    for (E_Int j = M->xcells[i]; j < M->xcells[i+1]; j++) {
      E_Int face = M->NFACE[j];
      if (M->owner[face] == -1) M->owner[face] = i;
      else M->neigh[face] = i;
    }
  }

  // TODO(Imad): replace with K_METRIC functions
  // Init mesh connectivity
  topo_init_mesh(M);

  // Parse comm data
  M->nppatches = PyList_Size(comm_data);
  M->ppatches = new proc_patch[M->nppatches];
  for (E_Int i = 0; i < M->nppatches; i++) {
    // Comm array
    PyObject *o = PyList_GetItem(comm_data, i);
    E_Int nfld, nneis;

    // Neighbour proc
    PyObject *nei_proc = PyList_GetItem(o, 0);
    M->ppatches[i].nei_proc = PyLong_AsLong(nei_proc);

    // Patch faces
    PyObject *farr = PyList_GetItem(o, 1);
    E_Int npfaces;
    K_NUMPY::getFromNumpyArray(farr, M->ppatches[i].faces, npfaces, nfld, false);
    assert(nfld == 1);
    M->ppatches[i].nfaces = npfaces;
    
    // Make them 0-based
    for (E_Int j = 0; j < npfaces; j++)
      M->ppatches[i].faces[j] -= 1;

    // Corresponding neighbours
    PyObject *narr = PyList_GetItem(o, 2);
    K_NUMPY::getFromNumpyArray(narr, M->ppatches[i].gneis, nneis, nfld, true);
    assert(nneis == M->ppatches[i].nfaces);
  }

  // Process solution fields
  E_Int csize = PyList_Size(solc);
  if (csize == 0) {
    RELEASESHAREDU(arr, f, cn);
    mesh_free(M);
    PyErr_SetString(PyExc_ValueError,
      "adaptMesh(): empty list of solution fields.");
    return NULL;
  }
  
  E_Int ok_sols = 1;
  E_Float **csols = (E_Float **)XCALLOC(csize, sizeof(E_Float *));
  for (E_Int i = 0; i < csize; i++) {
    PyObject *csol = PyList_GetItem(solc, i);
    E_Int nfld, nc;
    ret = K_NUMPY::getFromNumpyArray(csol, csols[i], nc, nfld, true);
    ok_sols = (ret == 1) && (nfld == 1) && (nc == M->ncells);
    if (!ok_sols) {
      PyErr_Format(PyExc_ValueError,
      "adaptMesh(): bad %d-th solution field.", i);
      break;
    }
  }

  if (!ok_sols) {
    mesh_free(M);
    for (E_Int i = 0; i < csize; i++)
      XFREE(csols[i]);
    XFREE(csols);
    RELEASESHAREDU(arr, f, cn);
    return NULL;
  }

  // Compute cell refinement levels
  ret = make_ref_data(M, csols, csize, freezeVector);
  if (ret == 1) {
    PyErr_SetString(PyExc_ValueError,
      "adaptMesh(): refinement cells not aligned with freezeVector. Aborting.");
    return NULL;
  }

  // Redistribute
  mesh *rM = redistribute_mesh(M);

  mesh_free(M);
  for (E_Int i = 0; i < csize; i++)
    XFREE(csols[i]);
  XFREE(csols);

  // adapt!
  if (rM->pid == 0)
    printf("Adapting...\n");
  
  MPI_Barrier(MPI_COMM_WORLD);
  clock_t tic = clock();

  // Isolate refinement cells
  std::vector<E_Int> ref_cells = get_ref_cells(rM);
  E_Int nref_cells = ref_cells.size();

  rM->ref_iter = 0;
  E_Int max_ref_iter = 10;

  tree ct(rM->ncells, 8);
  tree ft(rM->nfaces, 4);
  // Resize data structures (isotropic resizing, faster)
  resize_data_for_refinement(rM, &ct, &ft, nref_cells);

  while (nref_cells > 0) {
    if (++rM->ref_iter > max_ref_iter)
      break;

    // Refine cells
    for (E_Int i = 0; i < nref_cells; i++) {
      E_Int cell = ref_cells[i];
      E_Int *pr = &rM->ref_data[3*cell];
      if (pr[0] && pr[1] && pr[2]) {
        cut_cell_xyz(cell, rM, &ct, &ft);
      } else if (pr[0] && pr[1] && !pr[2]) {
        cut_cell_xy(cell, rM, &ct, &ft);
      } else if (pr[0] && !pr[1] && pr[2]) {
        cut_cell_xz(cell, rM, &ct, &ft);
      } else if (!pr[0] && pr[1] && pr[2]) {
        cut_cell_yz(cell, rM, &ct, &ft);
      } else if (pr[0] && !pr[1] && !pr[2]) {
        cut_cell_x(cell, rM, &ct, &ft);
      } else if (!pr[0] && pr[1] && !pr[2]) {
        cut_cell_y(cell, rM, &ct, &ft);
      } else if (!pr[0] && !pr[1] && pr[2]) {
        cut_cell_z(cell, rM, &ct, &ft);
      } else {
        assert(0);
      }
    }

    mesh_save_memory(rM);
    //tree_save_memory(&ct, rM->ncells);
    //tree_save_memory(&ft, rM->nfaces);

    // Update refinement data
    rM->ref_data = (E_Int *)XRESIZE(rM->ref_data, 3*rM->ncells * sizeof(E_Int));
    std::vector<E_Int> new_ref_cells(8*nref_cells);
    E_Int nref_cells_next_gen = 0;

    for (E_Int i = 0; i < nref_cells; i++) {
      E_Int cell = ref_cells[i];
      assert(ct.enabled[cell] == 0);
      E_Int *pr = &rM->ref_data[3*cell];
      for (E_Int j = 0; j < 3; j++)
        pr[j] = std::max((E_Int)0, pr[j]-1);

      if (is_cell_to_refine(pr)) {
        E_Int nchildren = tree_get_nchildren(&ct, cell);
        E_Int *children = tree_get_children(&ct, cell);
        for (E_Int j = 0; j < nchildren; j++) {
          assert(ct.enabled[children[j]] == 1);
          E_Int *prc = &rM->ref_data[3*children[j]];
          for (E_Int k = 0; k < 3; k++)
            prc[k] = pr[k];
          new_ref_cells[nref_cells_next_gen++] = children[j];
        }
      }

      for (E_Int j = 0; j < 3; j++)
        pr[j] = 0;
    }

    if (nref_cells_next_gen == 0) {
      break;
    }
    
    ref_cells.resize(nref_cells_next_gen);
    for (E_Int i = 0; i < nref_cells_next_gen; i++)
      ref_cells[i] = new_ref_cells[i];

    nref_cells = nref_cells_next_gen;

    resize_data_for_refinement(rM, &ct, &ft, nref_cells);
  }

  MPI_Barrier(MPI_COMM_WORLD);
  clock_t toc = clock();
  E_Float adapt_time = ((E_Float)(toc-tic)) / CLOCKS_PER_SEC;

  E_Int gnc = 0;
  MPI_Allreduce(&rM->ncells, &gnc, 1, XMPI_INT, MPI_SUM, MPI_COMM_WORLD);

  E_Float meshmem = mesh_memsize(rM);
  E_Float ftmem = tree_memsize(&ft);
  E_Float ctmem = tree_memsize(&ct);
  E_Float lmem = meshmem + ftmem + ctmem;
  E_Float gmem = 0;
  MPI_Allreduce(&lmem, &gmem, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  
  if (rM->pid == 0) {
    printf("Adaptation time: %.2f s\n", adapt_time);
    printf("Total number of cells: " SF_D_ "\n", gnc);
    printf("Memory: %.1f MB\n", gmem);
  }

  // output
  if (rM->pid == 0)
    puts("Exporting...");

  MPI_Barrier(MPI_COMM_WORLD);
  tic = clock();

  mesh_leaves cM = mesh_get_leaves(rM, &ct, &ft, conformize);

  const char *varStringOut = "CoordinateX,CoordinateY,CoordinateZ";

  PyObject* m = K_ARRAY::buildArray3(3, varStringOut, cM.nepoints, cM.necells,
    cM.nefaces, "NGON", cM.XFACES[cM.nefaces], cM.XCELLS[cM.necells], 3, 
    false, 3);

  FldArrayF *fo;
  FldArrayI *cno;
  K_ARRAY::getFromArray3(m, fo, cno);

  for (E_Int n = 0; n < 3; n++) {
    E_Float *pt = fo->begin(n+1);
    for (E_Int i = 0; i < cM.nepoints; i++)
      pt[i] = cM.XYZ[3*i+n];
  }

  E_Int *ngono = cno->getNGon();
  E_Int *nfaceo = cno->getNFace();
  E_Int *indPGo = cno->getIndPG();
  E_Int *indPHo = cno->getIndPH();

  for (E_Int i = 0; i <= cM.nefaces; i++) indPGo[i] = cM.XFACES[i];

  for (E_Int i = 0; i <= cM.necells; i++) indPHo[i] = cM.XCELLS[i];

  E_Int *ptr = ngono;
  E_Int start, end;

  for (E_Int i = 0; i < cM.nefaces; i++)
  {
    start = cM.XFACES[i];
    end = cM.XFACES[i+1];
    for (E_Int j = start; j < end; j++)
    { *ptr = cM.epoints[cM.NGON[j]]+1; ptr++; }
  }

  ptr = nfaceo;
  for (E_Int i = 0; i < cM.necells; i++)
  {
    start = cM.XCELLS[i];
    end = cM.XCELLS[i+1];
    for (E_Int j = start; j < end; j++)
    { *ptr = cM.efaces[cM.NFACE[j]]+1; ptr++; }
  }

  E_Int nleaves = 0;
  MPI_Allreduce(&cM.necells, &nleaves, 1, XMPI_INT, MPI_SUM, MPI_COMM_WORLD);

  toc = clock();
  E_Float extract_time = ((E_Float)(toc-tic)) / CLOCKS_PER_SEC;
  if (rM->pid == 0)
    printf("Extracted " SF_D_ " leaves in %.2f s\n", nleaves, extract_time);

  mesh_free(rM);

  return m;
}
