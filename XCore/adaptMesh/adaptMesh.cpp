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

PyObject *K_XCORE::adaptMesh(PyObject *self, PyObject *args)
{
  PyObject *arr, *comm_data, *solc, *gcells, *gfaces, *gpoints;

  if (!PyArg_ParseTuple(args, "OOOOOO", &arr, &comm_data, &solc, &gcells, &gfaces, &gpoints)) {
    PyErr_SetString(PyExc_ValueError, "adaptMesh(): wrong input.");
    return NULL;
  }

  FldArrayI *cn;
  FldArrayF *f;
  E_Int res = K_ARRAY::getFromArray3(arr, f, cn);

  // TODO(Imad): more input error checking
  if (res != 2) {
    PyErr_SetString(PyExc_TypeError, "adaptMesh(): input array should be NGON2.");
    return NULL;
  }

  // init mesh
  mesh *M = new mesh();
  M->NGON = cn->getNGon();
  M->NFACE = cn->getNFace();
  M->xfaces = cn->getIndPG();
  M->xcells = cn->getIndPH();

  // coordinates
  M->npoints = f->getSize();
  M->xyz = (E_Float *)malloc(3*M->npoints * sizeof(E_Float));
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

  // zero based data
  shift_data(M);
  
  // init parent elements
  M->owner = (E_Int *)malloc(M->nfaces * sizeof(E_Int));
  M->neigh = (E_Int *)malloc(M->nfaces * sizeof(E_Int));
  memset(M->owner, -1, M->nfaces*sizeof(E_Int));
  memset(M->neigh, -1, M->nfaces*sizeof(E_Int));
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

    E_Int *sorted_pfaces = (E_Int *)malloc(npfaces * sizeof(E_Int));
    E_Int *sorted_gneis = (E_Int *)malloc(npfaces * sizeof(E_Int));
    for (E_Int j = 0; j < npfaces; j++) {
      sorted_pfaces[j] = M->ppatches[i].faces[indices[j]];
      sorted_gneis[j] = M->ppatches[i].gneis[indices[j]];
    }

    free(M->ppatches[i].faces);
    free(M->ppatches[i].gneis);

    M->ppatches[i].faces = sorted_pfaces;
    M->ppatches[i].gneis = sorted_gneis;

    // replace with local ids
    for (E_Int j = 0; j < npfaces; j++) {
      assert(M->FT.find(M->ppatches[i].faces[j]) != M->FT.end());
      M->ppatches[i].faces[j] = M->FT[M->ppatches[i].faces[j]];
    }
  }

  // process solution fields
  E_Int csize = PyList_Size(solc);
  if (csize == 0) {
    PyErr_SetString(PyExc_ValueError, "adaptMesh(): empty solution.");
    return NULL;
  }
  
  E_Float **csols = (E_Float **)malloc(csize * sizeof(E_Float *));
  for (E_Int i = 0; i < csize; i++) {
    PyObject *csol = PyList_GetItem(solc, i);
    E_Int nfld;
    res = K_NUMPY::getFromNumpyArray(csol, csols[i], M->ncells, nfld, true); 
    assert(res == 1);
  }

  // TODO(Imad): for now, assume only one solution field. Later, implement metric interpolation
  
  // compute hessians
  E_Float *H = compute_hessian(M, csols[0]);
  
  // process hessians
  hessian_to_metric(H, M);

  // compute refinement data
  compute_ref_data(M, H);

  // process refinement data
  smooth_ref_data(M);

  // redistribute
  //mesh *rM = redistribute_mesh(M);

  // isolate refinement cells
  E_Int nref_cells = -1;
  E_Int nref_faces = -1;
  std::vector<E_Int> ref_cells = get_ref_cells(M, &nref_cells, &nref_faces);

  E_Int ref_iter = 0;
  E_Int max_ref_iter = 10;

  tree ct(M->ncells, 8);
  tree ft(M->nfaces, 4);

  while (nref_cells > 0) {
    if (++ref_iter > max_ref_iter)
      break;

    // resize data structures (isotropic resizing, faster)
    resize_data_for_refinement(M, &ct, &ft, nref_cells, nref_faces);

    // refine cells
    for (E_Int i = 0; i < nref_cells; i++) {
      E_Int cell = ref_cells[i];
      E_Int *pr = &M->ref_data[3*cell];
      if (pr[0] && pr[1] && pr[2]) {
        cut_cell_xyz(cell, M, &ct, &ft);
      } else if (pr[0] && pr[1] && !pr[2]) {
        cut_cell_xy(cell, M, &ct, &ft);
      } else if (pr[0] && !pr[1] && pr[2]) {
        cut_cell_xz(cell, M, &ct, &ft);
      } else if (!pr[0] && pr[1] && pr[2]) {
        cut_cell_yz(cell, M, &ct, &ft);
      } else if (pr[0] && !pr[1] && !pr[2]) {
        cut_cell_x(cell, M, &ct, &ft);
      } else if (!pr[0] && pr[1] && !pr[2]) {
        cut_cell_y(cell, M, &ct, &ft);
      } else if (!pr[0] && !pr[1] && pr[2]) {
        cut_cell_z(cell, M, &ct, &ft);
      } else {
        assert(0);
      }
    }

    // update refinement data
    M->ref_data = (E_Int *)realloc(M->ref_data, 3*M->ncells * sizeof(E_Int));
    std::vector<E_Int> new_ref_cells(8*nref_cells);
    E_Int nref_cells_next_gen = 0;

    for (E_Int i = 0; i < nref_cells; i++) {
      E_Int cell = ref_cells[i];
      E_Int *pr = &M->ref_data[3*cell];
      for (E_Int j = 0; j < 3; j++)
        pr[j] = std::max(0, pr[j]-1);

      if (is_cell_to_refine(pr)) {
        E_Int nchildren = tree_get_nchildren(&ct, cell);
        E_Int *children = tree_get_children(&ct, cell);
        for (E_Int j = 0; j < nchildren; j++) {
          E_Int *prc = &M->ref_data[3*children[j]];
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

    nref_faces = 6*nref_cells;
    
    ref_cells.resize(nref_cells_next_gen);
    for (E_Int i = 0; i < nref_cells_next_gen; i++)
      ref_cells[i] = new_ref_cells[i];

    nref_cells = nref_cells_next_gen;
  }

  //MPI_Barrier(MPI_COMM_WORLD);
  E_Int gnc;
  MPI_Allreduce(&M->ncells, &gnc, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  if (M->pid == 0)
    printf("Final number of cells: %d\n", gnc);

  free(H);

  // output
  const char *varString = "CoordinateX,CoordinateY,CoordinateZ";
  PyObject* m = K_ARRAY::buildArray3(3, varString, M->npoints, M->ncells, M->nfaces, "NGON",
                               M->xfaces[M->nfaces], M->xcells[M->ncells], 3,
                               false, 3);

  FldArrayF *fo;
  FldArrayI *cno;
  K_ARRAY::getFromArray3(m, fo, cno);

  assert(M->xfaces[M->nfaces] == 4*M->nfaces);
  assert(M->xcells[M->ncells] == 6*M->ncells);

  for (E_Int n = 0; n < 3; n++) {
    E_Float *pt = fo->begin(n+1);
    for (E_Int i = 0; i < M->npoints; i++)
      pt[i] = M->xyz[3*i+n];
  }

  E_Int * ngono = cno->getNGon();
  E_Int * nfaceo = cno->getNFace();
  E_Int * indPGo = cno->getIndPG();
  E_Int * indPHo = cno->getIndPH();
  for (E_Int i = 0; i <= M->nfaces; i++) indPGo[i] = M->xfaces[i];
  for (E_Int i = 0; i <= M->ncells; i++) indPHo[i] = M->xcells[i];
  E_Int* ptr = ngono;
  E_Int start, end;
  for (E_Int i = 0; i < M->nfaces; i++)
  {
    start = M->xfaces[i];
    end = M->xfaces[i+1];
    for (E_Int j = start; j < end; j++)
    { *ptr = M->NGON[j]+1; ptr++; }
  }
  ptr = nfaceo;
  for (E_Int i = 0; i < M->ncells; i++)
  {
    start = M->xcells[i];
    end = M->xcells[i+1];
    for (E_Int j = start; j < end; j++)
    { *ptr = M->NFACE[j]+1; ptr++; }
  }

  return m;
}
