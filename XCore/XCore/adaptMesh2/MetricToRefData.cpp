#include "Proto.h"
#include <stack>

#define MINREF -10

static
void make_ref_data_hexa(E_Int cell, AMesh *M, E_Float *pM, const pDirs &Dirs)
{
  E_Float L0, L1, L2, dd[3];
  //E_Float l[3];

  // Compute length in metric space

  K_MATH::sym3mat_dot_vec(pM, Dirs.I, dd);
  L0 = K_MATH::norm(dd, 3);

  K_MATH::sym3mat_dot_vec(pM, Dirs.J, dd);
  L1 = K_MATH::norm(dd, 3);

  K_MATH::sym3mat_dot_vec(pM, Dirs.K, dd);
  L2 = K_MATH::norm(dd, 3);

  M->ref_data[cell] = 0;

  if (L0 >= M->Tr || L1 >= M->Tr || L2 >= M->Tr) {
    M->ref_data[cell] = 1;
  }/* else if (L0 <= M->Tu && L1 <= M->Tu && L2 <= M->Tu) {
    if (M->cellTree->level(cell) != 0)
      M->ref_data[cell] = -1;
  }*/
}

void smooth_ref_data(AMesh *M)
{
    // Init refinement stack
    std::stack<E_Int> stk;
    for (E_Int i = 0; i < M->ncells; i++) {
      if (M->ref_data[i] > 0)
        stk.push(i);
    }

    // Smooth locally
    while (!stk.empty()) {
      E_Int cell = stk.top();
      stk.pop();

      E_Int nf = -1;
      E_Int pf[24];
      get_full_cell(cell, M, nf, pf);

      std::vector<E_Int> neis(nf);
      for (E_Int i = 0; i < nf; i++)
        neis[i] = get_neighbour(cell, pf[i], M);

      E_Int incr_cell = M->ref_data[cell] + M->cellTree->level(cell);

      for (E_Int i = 0; i < nf; i++) {
        E_Int nei = neis[i];
        if (nei == -1) continue;

        E_Int incr_nei = M->ref_data[nei] + M->cellTree->level(nei);

        E_Int diff = abs(incr_nei - incr_cell);

        if (diff <= 1) continue;

        E_Int cell_to_mod = incr_cell > incr_nei ? nei : cell;

        E_Int &rmod = M->ref_data[cell_to_mod];

        rmod += 1;
        
        stk.push(cell_to_mod);
      }
    }

}


void smooth_ref_data_parallel(AMesh *M)
{
  E_Int exchange, max_exchanges;
  exchange = 0;
  max_exchanges = 10;

  for (E_Int i = 0; i < M->npatches; i++) {
    Patch *P = &M->patches[i];
    P->sbuf_i = (E_Int *)XRESIZE(P->sbuf_i, P->nf*sizeof(E_Int));
    P->rbuf_i = (E_Int *)XRESIZE(P->rbuf_i, P->nf*sizeof(E_Int));
  }

  while (++exchange <= max_exchanges) {
    
    // Init refinement stack
    std::stack<E_Int> stk;
    for (E_Int i = 0; i < M->ncells; i++) {
      if (M->ref_data[i] > 0)
        stk.push(i);
    }

    // Smooth locally
    while (!stk.empty()) {
      E_Int cell = stk.top();
      stk.pop();

      E_Int nf = -1;
      E_Int pf[24];
      get_full_cell(cell, M, nf, pf);

      std::vector<E_Int> neis(nf);
      for (E_Int i = 0; i < nf; i++)
        neis[i] = get_neighbour(cell, pf[i], M);

      E_Int incr_cell = M->ref_data[cell] + M->cellTree->level(cell);

      for (E_Int i = 0; i < nf; i++) {
        E_Int nei = neis[i];
        if (nei == -1) continue;

        E_Int incr_nei = M->ref_data[nei] + M->cellTree->level(nei);

        E_Int diff = abs(incr_nei - incr_cell);

        if (diff <= 1) continue;

        E_Int cell_to_mod = incr_cell > incr_nei ? nei : cell;

        E_Int &rmod = M->ref_data[cell_to_mod];

        rmod += 1;
        
        stk.push(cell_to_mod);
      }
    }

    // Exchange proc ref data
    MPI_Barrier(MPI_COMM_WORLD);
    E_Int lstop = 0; // No more smoothing required locally

    for (E_Int i = 0; i < M->npatches; i++) {
      Patch *P = &M->patches[i];

      for (E_Int j = 0; j < P->nf; j++) {
        E_Int own = M->owner[P->pf[j]];
        P->sbuf_i[j] = M->ref_data[own] + M->cellTree->level(own);
      }

      MPI_Isend(P->sbuf_i, P->nf, XMPI_INT, P->nei, M->pid, MPI_COMM_WORLD, &M->req[M->nrq++]); 
      MPI_Irecv(P->rbuf_i, P->nf, XMPI_INT, P->nei, P->nei, MPI_COMM_WORLD, &M->req[M->nrq++]); 
    }
    Comm_waitall(M);


    for (E_Int i = 0; i < M->npatches; i++) {
      Patch *P = &M->patches[i];

      for (E_Int j = 0; j < P->nf; j++) {
        E_Int own = M->owner[P->pf[j]];
        E_Int oval = M->ref_data[own] + M->cellTree->level(own);
        E_Int nval = P->rbuf_i[j];

        //E_Int diff = abs(nval-oval);

        if (nval > oval + 1) {
          M->ref_data[own] += 1;
          lstop = 1;
        }
      }
    }
    
    E_Int gstop;
    MPI_Allreduce(&lstop, &gstop, 1, XMPI_INT, MPI_SUM, MPI_COMM_WORLD);

    if (gstop == 0) // No more smoothing required globally
      break;
  }

  //printf("%d -> exchanged %d times\n", M->pid, exchange);

  if (exchange > max_exchanges)
    fprintf(stderr, "Warning: smoothing exceeded max_exchanges!\n");
}

static
void make_ref_data_tetra(E_Int cell, AMesh *M, E_Float *pM, const pDirs &Dirs)
{}

static
void make_ref_data_penta(E_Int cell, AMesh *M, E_Float *pM, const pDirs &Dirs)
{}

static
void make_ref_data_pyra(E_Int cell, AMesh *M, E_Float *pM, const pDirs &Dirs)
{}

static
void compute_ref_data(AMesh *M, E_Float *metric, const std::vector<pDirs> &Dirs)
{
  for (E_Int i = 0; i < M->ncells; i++) {
    E_Int type = M->cellTree->type(i);
    E_Float *pM = &metric[6*i];
    switch (type) {
      case HEXA:
        make_ref_data_hexa(i, M, pM, Dirs[i]);
        break;
      case TETRA:
        make_ref_data_tetra(i, M, pM, Dirs[i]);
        break;
      case PENTA:
        make_ref_data_penta(i, M, pM, Dirs[i]);
        break;
      case PYRA:
        make_ref_data_pyra(i, M, pM, Dirs[i]);
        break;
      default:
        assert(0);
        break;
    }
  }

  if (!M->unrefine) {
    for (E_Int i = 0; i < M->ncells; i++) {
      if (M->ref_data[i] < 0)
        M->ref_data[i] = 0;
    }
  }
}

static
void make_pdirs_tetra(E_Int cell, AMesh *M, pDirs &Dirs)
{}

static
void make_pdirs_penta(E_Int cell, AMesh *M, pDirs &Dirs)
{}

static
void make_pdirs_pyra(E_Int cell, AMesh *M, pDirs &Dirs)
{}

void reconstruct_parent_quad(E_Int face, AMesh *M, E_Int pn[4])
{
  Children *children = M->faceTree->children(face);
  if (children == NULL) {
    memcpy(pn, &M->ngon[M->indPG[face]], 4*sizeof(E_Int));
  } else {
    if (children->n == 2) {
      E_Int *pn0 = get_facets(children->pc[0], M->ngon, M->indPG);
      E_Int *pn1 = get_facets(children->pc[1], M->ngon, M->indPG);
      E_Int dir = pn0[1] == pn1[0] ? DIRX : DIRY;
      if (dir == DIRX) {
        pn[0] = pn0[0]; pn[1] = pn1[1]; pn[2] = pn1[2]; pn[3] = pn0[3];
      } else {
        pn[0] = pn0[0]; pn[1] = pn0[1]; pn[2] = pn1[2]; pn[3] = pn1[3];
      }
    } else {
      for (E_Int i = 0; i < children->n; i++)
        pn[i] = get_facets(children->pc[i], M->ngon, M->indPG)[i];
    }
  }
}

static
void make_pdirs_hexa(E_Int cell, AMesh *M, pDirs &Dirs)
{
  E_Int *pf = &M->nface[M->indPH[cell]];

  E_Int p0[4], p1[4];
  E_Float f0[3], f1[3];

  // LFT && RGT
  reconstruct_parent_quad(pf[2], M, p0);
  reconstruct_parent_quad(pf[3], M, p1);
  f0[0] = f0[1] = f0[2] = 0.0;
  f1[0] = f1[1] = f1[2] = 0.0;

  for (E_Int i = 0; i < 4; i++) {
    f0[0] += M->x[p0[i]];
    f0[1] += M->y[p0[i]];
    f0[2] += M->z[p0[i]];
    f1[0] += M->x[p1[i]];
    f1[1] += M->y[p1[i]];
    f1[2] += M->z[p1[i]];
  }

  for (E_Int i = 0; i < 3; i++) Dirs.I[i] = 0.25*(f1[i] - f0[i]);

  // FRO && BCK
  reconstruct_parent_quad(pf[4], M, p0);
  reconstruct_parent_quad(pf[5], M, p1);
  f0[0] = f0[1] = f0[2] = 0.0;
  f1[0] = f1[1] = f1[2] = 0.0;

  for (E_Int i = 0; i < 4; i++) {
    f0[0] += M->x[p0[i]];
    f0[1] += M->y[p0[i]];
    f0[2] += M->z[p0[i]];
    f1[0] += M->x[p1[i]];
    f1[1] += M->y[p1[i]];
    f1[2] += M->z[p1[i]];
  }

  for (E_Int i = 0; i < 3; i++) Dirs.J[i] = 0.25*(f1[i] - f0[i]);

  // LFT && RGT
  reconstruct_parent_quad(pf[0], M, p0);
  reconstruct_parent_quad(pf[1], M, p1);
  f0[0] = f0[1] = f0[2] = 0.0;
  f1[0] = f1[1] = f1[2] = 0.0;

  for (E_Int i = 0; i < 4; i++) {
    f0[0] += M->x[p0[i]];
    f0[1] += M->y[p0[i]];
    f0[2] += M->z[p0[i]];
    f1[0] += M->x[p1[i]];
    f1[1] += M->y[p1[i]];
    f1[2] += M->z[p1[i]];
  }

  for (E_Int i = 0; i < 3; i++) Dirs.K[i] = 0.25*(f1[i] - f0[i]);
}


static
void compute_principal_vecs(AMesh *M, std::vector<pDirs> &Dirs)
{
  for (E_Int i = 0; i < M->ncells; i++) {
    switch (M->cellTree->type(i)) {
      case TETRA:
        make_pdirs_tetra(i, M, Dirs[i]);
        break;
      case PENTA:
        make_pdirs_penta(i, M, Dirs[i]);
        break;
      case PYRA:
        make_pdirs_pyra(i, M, Dirs[i]);
        break;
      case HEXA:
        make_pdirs_hexa(i, M, Dirs[i]);
        break;
      default:
        assert(0);
        break;
    }
  }
}

PyObject *K_XCORE::_metricToRefData(PyObject *self, PyObject *args)
{
  PyObject *METRIC, *AMESH;

  if (!PYPARSETUPLE_(args, OO_, &METRIC, &AMESH)) {
    RAISE("Wrong input.");
    return NULL;
  }

  AMesh *M = (AMesh *)PyCapsule_GetPointer(AMESH, "AMesh");

  // Parse field
  E_Float *metric = NULL;
  E_Int nfld, size;
  K_NUMPY::getFromNumpyArray(METRIC, metric, size, nfld, true);
  //assert(ret == 1 && nfld == 1 && size == M->ncells*6);

  // Init ref data
  M->ref_data = (int *)XRESIZE(M->ref_data, M->ncells*sizeof(int));
  for (E_Int i = 0; i < M->ncells; i++) M->ref_data[i] = 0;

  // Compute cells principal directions
  std::vector<pDirs> Dirs(M->ncells);
  compute_principal_vecs(M, Dirs);

  // Compute ref data
  compute_ref_data(M, metric, Dirs);

  // Smooth ref data
  smooth_ref_data_parallel(M);

  return Py_None;
}

static
E_Float hessian_norm(E_Float M[6])
{
  E_Float row0 = fabs(M[0]) + fabs(M[1]) + fabs(M[2]);
  E_Float row1 = fabs(M[1]) + fabs(M[3]) + fabs(M[4]);
  E_Float row2 = fabs(M[2]) + fabs(M[4]) + fabs(M[5]);

  return std::max(row0, std::max(row1, row2));
}

static
E_Float gradient_norm(E_Float G[3])
{
  return fabs(G[0]) + fabs(G[1]) + fabs(G[2]);
}

PyObject *K_XCORE::_makeRefDataFromGradAndHess(PyObject *self, PyObject *args)
{
  PyObject *AMESH, *FIELD, *GRAD, *HESS;

  if (!PYPARSETUPLE_(args, OOOO_, &AMESH, &FIELD, &GRAD, &HESS)) {
    RAISE("Wrong input.");
    return NULL;
  }

  AMesh *M = (AMesh *)PyCapsule_GetPointer(AMESH, "AMesh");

  // Parse field/grad/hess
  E_Float *field, *grad, *hess;
  E_Int nfld, size, ret;
  ret = K_NUMPY::getFromNumpyArray(FIELD, field, size, nfld, true);
  assert(ret == 1 && nfld == 1 && size == M->ncells);
  ret = K_NUMPY::getFromNumpyArray(GRAD, grad, size, nfld, true);
  assert(ret == 1 && nfld == 1 && size == M->ncells*3);
  ret = K_NUMPY::getFromNumpyArray(HESS, hess, size, nfld, true);
  assert(ret == 1 && nfld == 1 && size == M->ncells*6);

  npy_intp dims[2];
  dims[0] = (npy_intp)M->ncells;
  dims[1] = 1;

  PyArrayObject *R = (PyArrayObject *)PyArray_SimpleNew(1, dims, E_NPY_INT);
  E_Int *ptr = (E_Int *)PyArray_DATA(R);
  memset(ptr, 0, M->ncells*sizeof(E_Int));

  // Compute cells principal directions
  std::vector<pDirs> Dirs(M->ncells);
  compute_principal_vecs(M, Dirs);

  // Compute ref data
  E_Float mean = 0.0;
  for (E_Int i = 0; i < M->ncells; i++) {
    mean += field[i];
  }

  E_Float gmean;
  MPI_Allreduce(&mean, &gmean, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  E_Int gncells;
  MPI_Allreduce(&M->ncells, &gncells, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

  mean = fabs(gmean / gncells);

  for (E_Int i = 0; i < M->ncells; i++) {
    
    // Ref length
    E_Float hi = K_MATH::norm(Dirs[i].I, 3);
    E_Float hj = K_MATH::norm(Dirs[i].J, 3);
    E_Float hk = K_MATH::norm(Dirs[i].K, 3);
    E_Float h = std::min(hi, std::min(hj, hk));

    E_Float val = h*h*hessian_norm(&hess[6*i]);
    val /= h*gradient_norm(&grad[3*i]) + M->eps*mean;

    if (val >= M->Tr && h > M->hmin) ptr[i] = 1;
  }

  return (PyObject *)R;
}
