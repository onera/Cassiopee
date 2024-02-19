#include "Proto.h"

E_Float gradient_norm_inf(E_Float G[3])
{
  return fmax(fabs(G[0]), fmax(fabs(G[1]), fabs(G[2])));
}

static
void make_A_grad_matrices(AMesh *M, E_Int *indPH, E_Int *count_neis,
  E_Int *owner, E_Int *neigh, E_Float *cx, E_Float *cy, E_Float *cz,
  E_Float *lsqG)
{
  E_Int ncells = M->ncells;
  E_Int nfaces = M->nfaces;
  memset(count_neis, 0, ncells*sizeof(E_Int));

  // Internal faces
  for (E_Int i = 0; i < nfaces; i++) {
    E_Int nei = neigh[i];
    if (nei == -1) continue;

    E_Int own = owner[i];
    
    E_Float d[3] = {cx[nei]-cx[own],
                    cy[nei]-cy[own],
                    cz[nei]-cz[own]};

    E_Float *lo = &lsqG[3*(indPH[own] + count_neis[own]++)];
    E_Float *ln = &lsqG[3*(indPH[nei] + count_neis[nei]++)];

    for (E_Int j = 0; j < 3; j++) {
      lo[j] = d[j];
      ln[j] = -d[j];
    }
  }

  for (E_Int i = 0; i < M->npatches; i++) {
    Patch *P = &M->patches[i];

    E_Float *it = P->rbuf_f;

    for (E_Int j = 0; j < P->nf; j++) {
      E_Int own = owner[P->pf[j]];
      
      E_Float d[3];
      d[0] = *it++ - cx[own];
      d[1] = *it++ - cy[own];
      d[2] = *it++ - cz[own];
      it++; // skip field

      E_Float *lo = &lsqG[3*(indPH[own] + count_neis[own]++)];

      for (E_Int k = 0; k < 3; k++)
        lo[k] = d[k];
    }
  }

  // TODO(Imad): boundary contributions
}

static
void deduce_lsq_grad_matrices(AMesh *M, E_Int *indPH, E_Int *count_neis,
  E_Float *lsqG, E_Float *lsqGG)
{
  E_Int ncells = M->ncells;

  for (E_Int i = 0; i < ncells; i++) {
    E_Float *pGG = &lsqGG[9*i];
    E_Float *pG = &lsqG[3*indPH[i]];
    E_Int idx = 0;
    for (E_Int j = 0; j < 3; j++) {
      for (E_Int k = 0; k < 3; k++) {
        pGG[idx] = 0.0;
        for (E_Int l = 0; l < count_neis[i]; l++)
          pGG[idx] += pG[l*3+j] * pG[l*3+k];
        idx++;
      }
    }
  }
}

static
void make_grad_RHS_vector(AMesh *M, E_Int *indPH, E_Int *count_neis,
  E_Float *field, E_Float *b, E_Int *owner, E_Int *neigh)
{
  E_Int ncells = M->ncells;
  memset(count_neis, 0, ncells*sizeof(E_Int));

  E_Int nfaces = M->nfaces;

  // Construct b vector
  for (E_Int i = 0; i < nfaces; i++) {
    E_Int nei = neigh[i];
    if (nei == -1) continue;

    E_Int own = owner[i];
    b[indPH[own] + count_neis[own]++] = field[nei] - field[own];
    b[indPH[nei] + count_neis[nei]++] = field[own] - field[nei];
  }

  for (E_Int i = 0; i < M->npatches; i++) {
    Patch *P = &M->patches[i];

    E_Float *it = P->rbuf_f;

    for (E_Int j = 0; j < P->nf; j++) {
      E_Int own = owner[P->pf[j]];
      
      it += 3; // skip coordinates
      
      b[indPH[own] + count_neis[own]++] = *it++ - field[own];
    }
  }

  // TODO(Imad): boundary faces contributions
}

static
void make_gradient(AMesh *M, E_Int *indPH, E_Int *count_neis, E_Float *b,
  E_Float *G, E_Float *lsqG, E_Float *lsqGG)
{
  E_Int ncells = M->ncells;
  E_Float B[3];

  for (E_Int i = 0; i < ncells; i++) {
    // Make B vector = tA.b
    E_Float *pG = &lsqG[3*indPH[i]];
    E_Float *pb = &b[indPH[i]];
    for (E_Int j = 0; j < 3; j++) {
      B[j] = 0.0;
      for (E_Int k = 0; k < count_neis[i]; k++)
        B[j] += pb[k] * pG[k*3 + j];
    }

    // Solve
    E_Int converged = K_LINEAR::BiCGStab(&lsqGG[9*i], B, 3, &G[3*i]);

    if (!converged) {
      fprintf(stderr, "make_gradient: cell %d not converged.\n", i);
      exit(1);
    }
  }
}

PyObject *K_XCORE::computeGradient(PyObject *self, PyObject *args)
{
  PyObject *AMESH, *FIELD, *CX, *CY, *CZ, *OWN, *NEI;
  if (!PYPARSETUPLE_(args, OOO_ OOO_ O_, &AMESH, &FIELD, &CX, &CY,
      &CZ, &OWN, &NEI)) {
    RAISE("Wrong input.");
    return NULL;
  }

  if (!PyCapsule_IsValid(AMESH, "AMesh")) {
    RAISE("Bad mesh capsule.");
    return NULL;
  }

  AMesh *M = (AMesh *)PyCapsule_GetPointer(AMESH, "AMesh");

  // Parse field
  E_Float *field = NULL;
  E_Int nfld, size, ret;
  E_Int ncells = M->ncells;
  ret = K_NUMPY::getFromNumpyArray(FIELD, field, size, nfld, true);
  assert(ret == 1 && nfld == 1 && size == ncells);

  // Parse cell centers
  E_Float *cx, *cy, *cz;
  K_NUMPY::getFromNumpyArray(CX, cx, size, nfld, true);
  K_NUMPY::getFromNumpyArray(CY, cy, size, nfld, true);
  K_NUMPY::getFromNumpyArray(CZ, cz, size, nfld, true);

  // Parse owner and neigh
  E_Int *owner, *neigh;
  K_NUMPY::getFromNumpyArray(OWN, owner, size, nfld, true);
  K_NUMPY::getFromNumpyArray(NEI, neigh, size, nfld, true);

  // Compute gradient

  E_Int sizeNFace, *indPH;
  if (M->closed_indPH) {
    sizeNFace = M->closed_indPH[M->ncells];
    indPH = M->closed_indPH;
  } else {
    sizeNFace = M->indPH[M->ncells];
    indPH = M->indPH;
  }

  E_Int *count_neis = (E_Int *)XMALLOC(ncells * sizeof(E_Int));

  E_Float *lsqG = (E_Float *)XMALLOC(3*sizeNFace * sizeof(E_Float));
  E_Float *lsqGG = (E_Float *)XMALLOC(9*ncells * sizeof(E_Float));

  // Pre-exchange
  for (E_Int i = 0; i < M->npatches; i++) {
    Patch *P = &M->patches[i];

    P->sbuf_f = (E_Float *)XRESIZE(P->sbuf_f, 4*P->nf*sizeof(E_Float));
    P->rbuf_f = (E_Float *)XRESIZE(P->rbuf_f, 4*P->nf*sizeof(E_Float));
  }

  for (E_Int i = 0; i < M->npatches; i++) {
    Patch *P = &M->patches[i];

    E_Float *it = P->sbuf_f;

    for (E_Int j = 0; j < P->nf; j++) {
      E_Int own = owner[P->pf[j]];
      *it++ = cx[own]; 
      *it++ = cy[own]; 
      *it++ = cz[own]; 
      *it++ = field[own];
    }

    MPI_Isend(P->sbuf_f, 4*P->nf, MPI_DOUBLE, P->nei, M->pid, MPI_COMM_WORLD, &M->req[M->nrq++]);
    MPI_Irecv(P->rbuf_f, 4*P->nf, MPI_DOUBLE, P->nei, P->nei, MPI_COMM_WORLD, &M->req[M->nrq++]);
  }
  Comm_waitall(M);

  make_A_grad_matrices(M, indPH, count_neis, owner, neigh, cx, cy, cz, lsqG);
  
  deduce_lsq_grad_matrices(M, indPH, count_neis, lsqG, lsqGG);

  E_Float *b = (E_Float *)XMALLOC(sizeNFace * sizeof(E_Float));

  make_grad_RHS_vector(M, indPH, count_neis, field, b, owner, neigh);

  npy_intp dims[2];
  dims[0] = (npy_intp)(ncells*3);
  dims[1] = 1;

  PyArrayObject *G = (PyArrayObject *)PyArray_SimpleNew(1, dims, NPY_DOUBLE);
  E_Float *pG = (E_Float *)PyArray_DATA(G);

  make_gradient(M, indPH, count_neis, b, pG, lsqG, lsqGG);

  /*
  for (E_Int i = 0; i < M->ncells; i++) {
    E_Float *ptr = &pG[3*i];
    for (E_Int j = 0; j < 3; j++)
      printf("%.2f ", ptr[j]);
    puts("");
  }
  */

  XFREE(b);
  XFREE(count_neis);
  XFREE(lsqG);
  XFREE(lsqGG);

  return (PyObject *)G;
}


