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

#define PARSE_DICT(dict, key) parse_dictionary((dict), (key), __FUNCTION__)

static
PyObject *parse_dictionary(PyObject *dict, const char *key, const char *func)
{
  PyObject *strObj;
#if PY_VERSION_HEX >= 0x03000000
  strObj = PyUnicode_FromString(key);
#else
  strObj = PyString_FromString(key);
#endif

  PyObject *obj = PyDict_GetItem(dict, strObj);
  if (obj == NULL) {
    PyErr_Format(PyExc_KeyError, "%s: Missing key %s in dictionary", func, key);
  }
  Py_DECREF(strObj);
  return obj;
}

PyObject *K_XCORE::adaptMesh(PyObject *self, PyObject *args)
{
  PyObject *arr, *comm_data, *solc, *gcells, *gfaces, *gpoints;
  E_Int Gmax, iso_mode, conformize;
  E_Float Tr;
  PyObject *adaptDict;

  if (!PyArg_ParseTuple(args, "OOOOOOO", &arr, &comm_data, &solc, &gcells, &gfaces, &gpoints, &adaptDict)) {
    PyErr_SetString(PyExc_ValueError, "adaptMesh(): wrong input.");
    return NULL;
  }

  assert(PyDict_Check(adaptDict));

  // parse dictionary
  PyObject *obj;

  // Gmax
  obj = PARSE_DICT(adaptDict, "Gmax");
  if (obj == NULL) return NULL;
  Gmax = PyLong_AsLong(obj);

  // Tr
  obj = PARSE_DICT(adaptDict, "Tr");
  if (obj == NULL) return NULL;
  Tr = PyFloat_AsDouble(obj);

  // iso_mode
  obj = PARSE_DICT(adaptDict, "iso_mode");
  if (obj == NULL) return NULL;
  iso_mode = PyLong_AsLong(obj);

  // conformize
  obj = PARSE_DICT(adaptDict, "conformize");
  if (obj == NULL) return NULL;
  conformize = PyLong_AsLong(obj);

  FldArrayI *cn;
  FldArrayF *f;
  E_Int res = K_ARRAY::getFromArray3(arr, f, cn);

  // TODO(Imad): more input error checking
  if (res != 2) {
    PyErr_SetString(PyExc_TypeError,
      "adaptMesh(): input array should be NGON2.");
    return NULL;
  }

  // init mesh
  mesh *M = new mesh();
  M->NGON = cn->getNGon();
  M->NFACE = cn->getNFace();
  M->xfaces = cn->getIndPG();
  M->xcells = cn->getIndPH();
  M->Gmax = Gmax;
  M->Tr = Tr;
  M->iso_mode = iso_mode;

  // coordinates
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

  // parse my global cells, faces and points
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
  MPI_Allreduce(&M->ncells, &gncells, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

  if (M->pid == 0)
    printf("Total number of cells: %d\n", gncells);

  shift_data(M);

  // init parent elements
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

  // init mesh connectivity
  topo_init_mesh(M);

  // parse comm data
  M->nppatches = PyList_Size(comm_data);
  //M->ppatches = (proc_patch *)XMALLOC(M->nppatches * sizeof(proc_patch));
  M->ppatches = new proc_patch[M->nppatches];
  for (E_Int i = 0; i < M->nppatches; i++) {
    // comm array
    PyObject *o = PyList_GetItem(comm_data, i);
    E_Int nfld, nneis;

    // neighbour proc
    PyObject *nei_proc = PyList_GetItem(o, 0);
    M->ppatches[i].nei_proc = PyLong_AsLong(nei_proc);

    // patch faces
    PyObject *farr = PyList_GetItem(o, 1);
    E_Int npfaces;
    K_NUMPY::getFromNumpyArray(farr, M->ppatches[i].faces, npfaces, nfld, true);
    M->ppatches[i].nfaces = npfaces;
   
    // corresponding neighbours
    PyObject *narr = PyList_GetItem(o, 2);
    K_NUMPY::getFromNumpyArray(narr, M->ppatches[i].gneis, nneis, nfld, true);
    assert(nneis == M->ppatches[i].nfaces);

    // indices
    std::vector<E_Int> indices(M->ppatches[i].nfaces);
    for (E_Int j = 0; j < npfaces; j++)
      indices[j] = j;

    std::sort(indices.begin(), indices.end(),
      [&](E_Int a, E_Int b)
      {
        return M->ppatches[i].faces[a] < M->ppatches[i].faces[b];
      });

    E_Int *sorted_pfaces = (E_Int *)XCALLOC(npfaces, sizeof(E_Int));
    E_Int *sorted_gneis = (E_Int *)XCALLOC(npfaces, sizeof(E_Int));
    for (E_Int j = 0; j < npfaces; j++) {
      sorted_pfaces[j] = M->ppatches[i].faces[indices[j]];
      sorted_gneis[j] = M->ppatches[i].gneis[indices[j]];
    }

    XFREE(M->ppatches[i].faces);
    XFREE(M->ppatches[i].gneis);

    M->ppatches[i].faces = sorted_pfaces;
    M->ppatches[i].gneis = sorted_gneis;

    // replace with local ids
    for (E_Int j = 0; j < npfaces; j++) {
      assert(M->FT.find(M->ppatches[i].faces[j]) != M->FT.end());
      assert(M->ppatches[i].faces[j] > 0);
      M->ppatches[i].faces[j] = M->FT[M->ppatches[i].faces[j]];
    }
  }

  // process solution fields
  E_Int csize = PyList_Size(solc);
  if (csize == 0) {
    PyErr_SetString(PyExc_ValueError, "adaptMesh(): empty solution.");
    return NULL;
  }
  
  E_Float **csols = (E_Float **)XCALLOC(csize, sizeof(E_Float *));
  for (E_Int i = 0; i < csize; i++) {
    PyObject *csol = PyList_GetItem(solc, i);
    E_Int nfld;
    res = K_NUMPY::getFromNumpyArray(csol, csols[i], M->ncells, nfld, true); 
    assert(res == 1);
  }

  // TODO(Imad): for now, assume only one solution field.
  // Later, implement metric interpolation
  
  // compute hessians
  E_Float *H = compute_hessian(M, csols[0]);
    
  // process hessians
  hessian_to_metric(H, M);
 
  // compute refinement data
  compute_ref_data(M, H);

  // process refinement data
  smooth_ref_data(M);

  // redistribute
  mesh *rM = redistribute_mesh(M);

  mesh_free(M);
  for (E_Int i = 0; i < csize; i++)
    XFREE(csols[i]);
  XFREE(csols);
  XFREE(H);

  // adapt!
  if (rM->pid == 0)
    printf("Adapting...\n");
  
  MPI_Barrier(MPI_COMM_WORLD);
  clock_t tic = clock();

  // isolate refinement cells
  std::vector<E_Int> ref_cells = get_ref_cells(rM);
  E_Int nref_cells = ref_cells.size();

  rM->ref_iter = 0;
  E_Int max_ref_iter = 10;

  tree ct(rM->ncells, 8);
  tree ft(rM->nfaces, 4);
  // resize data structures (isotropic resizing, faster)
  resize_data_for_refinement(rM, &ct, &ft, nref_cells);

  while (nref_cells > 0) {
    if (++rM->ref_iter > max_ref_iter)
      break;

    // refine cells
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

    // update refinement data
    rM->ref_data = (E_Int *)XRESIZE(rM->ref_data, 3*rM->ncells * sizeof(E_Int));
    std::vector<E_Int> new_ref_cells(8*nref_cells);
    E_Int nref_cells_next_gen = 0;

    for (E_Int i = 0; i < nref_cells; i++) {
      E_Int cell = ref_cells[i];
      assert(ct.enabled[cell] == 0);
      E_Int *pr = &rM->ref_data[3*cell];
      for (E_Int j = 0; j < 3; j++)
        pr[j] = std::max(0, pr[j]-1);

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
  MPI_Allreduce(&rM->ncells, &gnc, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

  E_Float meshmem = mesh_memsize(rM);
  E_Float ftmem = tree_memsize(&ft);
  E_Float ctmem = tree_memsize(&ct);
  E_Float lmem = meshmem + ftmem + ctmem;
  E_Float gmem = 0;
  MPI_Allreduce(&lmem, &gmem, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  
  if (rM->pid == 0) {
    printf("Adaptation time: %.2f s\n", adapt_time);
    printf("Total number of cells: %d\n", gnc);
    printf("Memory: %.1f MB\n", gmem);
  }

  // output
  if (rM->pid == 0)
    puts("Exporting...");

  MPI_Barrier(MPI_COMM_WORLD);
  tic = clock();

  mesh_leaves cM = mesh_get_leaves(rM, &ct, &ft, conformize);

  const char *varString = "CoordinateX,CoordinateY,CoordinateZ";

  PyObject* m = K_ARRAY::buildArray3(3, varString, cM.nepoints, cM.necells,
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
  MPI_Allreduce(&cM.necells, &nleaves, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

  toc = clock();
  E_Float extract_time = ((E_Float)(toc-tic)) / CLOCKS_PER_SEC;
  if (rM->pid == 0)
    printf("Extracted %d leaves in %.2f s\n", nleaves, extract_time);

  mesh_free(rM);

  return m;
}
