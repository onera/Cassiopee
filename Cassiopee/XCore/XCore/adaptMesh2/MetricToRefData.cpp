#include "Proto.h"
#include <stack>

#define MINREF -10

static
void smooth_ref_data_local(AMesh *M)
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

void smooth_ref_data(AMesh *M)
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
    
    smooth_ref_data_local(M);

    // Exchange proc ref data
    MPI_Barrier(MPI_COMM_WORLD);
    E_Int lstop = 0; // No more smoothing required locally

    for (E_Int i = 0; i < M->npatches; i++) {
      Patch *P = &M->patches[i];

      for (E_Int j = 0; j < P->nf; j++) {
        E_Int own = M->owner[P->pf[j]];
        P->sbuf_i[j] = M->ref_data[own] + M->cellTree->level(own);
      }

      MPI_Isend(P->sbuf_i, P->nf, XMPI_INT, P->nei, M->pid, MPI_COMM_WORLD,
        &M->req[M->nrq++]); 
      MPI_Irecv(P->rbuf_i, P->nf, XMPI_INT, P->nei, P->nei, MPI_COMM_WORLD,
        &M->req[M->nrq++]); 
    }

    Comm_waitall(M);

    for (E_Int i = 0; i < M->npatches; i++) {
      Patch *P = &M->patches[i];

      for (E_Int j = 0; j < P->nf; j++) {
        E_Int own = M->owner[P->pf[j]];
        E_Int oval = M->ref_data[own] + M->cellTree->level(own);
        E_Int nval = P->rbuf_i[j];

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

  if (exchange > max_exchanges)
    fprintf(stderr, "Warning: smoothing exceeded max_exchanges!\n");
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

    E_Float val = h*h*hessian_norm_inf(&hess[6*i]);
    val /= h*gradient_norm_inf(&grad[3*i]) + M->eps*mean;

    if (val >= M->Tr && h > M->hmin) ptr[i] = 1;
  }

  return (PyObject *)R;
}