#include "Proto.h"

static
void make_A_grad_matrices(K_FLD::FldArrayI &cn, E_Int *count_neis,
  E_Int *owner, E_Int *neigh, E_Float *cx, E_Float *cy, E_Float *cz,
  E_Float *lsqG)
{
  E_Int ncells = cn.getNElts();
  E_Int nfaces = cn.getNFaces();
  memset(count_neis, 0, ncells*sizeof(E_Int));

  E_Int *indPH = cn.getIndPH();

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

  // TODO(Imad): boundary contributions
  // TODO(Imad): patch contributions
}

static
void deduce_lsq_grad_matrices(K_FLD::FldArrayI &cn, E_Int *count_neis,
  E_Float *lsqG, E_Float *lsqGG)
{
  E_Int ncells = cn.getNElts();
  E_Int *indPH = cn.getIndPH();

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
void make_grad_RHS_vector(K_FLD::FldArrayI &cn, E_Int *count_neis,
  E_Float *field, E_Float *b, E_Int *owner, E_Int *neigh)
{
  E_Int ncells = cn.getNElts();
  memset(count_neis, 0, ncells*sizeof(E_Int));

  E_Int nfaces = cn.getNFaces();
  E_Int *indPH = cn.getIndPH();

  // Construct b vector
  for (E_Int i = 0; i < nfaces; i++) {
    E_Int nei = neigh[i];
    if (nei == -1) continue;

    E_Int own = owner[i];
    b[indPH[own] + count_neis[own]++] = field[nei] - field[own];
    b[indPH[nei] + count_neis[nei]++] = field[own] - field[nei];
  }

  // TODO(Imad): boundary faces contributions
  // TODO(Imad): proc contributions
}

static
void make_gradient(K_FLD::FldArrayI &cn, E_Int *count_neis, E_Float *b,
  E_Float *G, E_Float *lsqG, E_Float *lsqGG)
{
  E_Int ncells = cn.getNElts();
  E_Int *indPH = cn.getIndPH();
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
  PyObject *ARRAY, *FIELD, *CX, *CY, *CZ, *OWN, *NEI;
  if (!PYPARSETUPLE_(args, OOO_ OOO_ O_, &ARRAY, &FIELD, &CX, &CY,
      &CZ, &OWN, &NEI)) {
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

  // Parse field
  E_Float *field = NULL;
  E_Int nfld, size;
  E_Int ncells = cn->getNElts();
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

  E_Int sizeNFace = cn->getSizeNFace();

  E_Int *count_neis = (E_Int *)XMALLOC(ncells * sizeof(E_Int));

  E_Float *lsqG = (E_Float *)XMALLOC(3*sizeNFace * sizeof(E_Float));
  E_Float *lsqGG = (E_Float *)XMALLOC(9*ncells * sizeof(E_Float));

  make_A_grad_matrices(*cn, count_neis, owner, neigh, cx, cy, cz, lsqG);
  
  deduce_lsq_grad_matrices(*cn, count_neis, lsqG, lsqGG);

  E_Float *b = (E_Float *)XMALLOC(sizeNFace * sizeof(E_Float));

  make_grad_RHS_vector(*cn, count_neis, field, b, owner, neigh);

  npy_intp dims[2];
  dims[0] = (npy_intp)(ncells*3);
  dims[1] = 1;

  PyArrayObject *G = (PyArrayObject *)PyArray_SimpleNew(1, dims, NPY_DOUBLE);
  E_Float *pG = (E_Float *)PyArray_DATA(G);

  make_gradient(*cn, count_neis, b, pG, lsqG, lsqGG);

  XFREE(b);
  XFREE(count_neis);
  XFREE(lsqG);
  XFREE(lsqGG);

  return (PyObject *)G;
}


