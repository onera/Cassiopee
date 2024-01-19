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

static
E_Float *compute_gradient(K_FLD::FldArrayI &cn, E_Float *field, E_Int *owner,
  E_Int *neigh, E_Float *cx, E_Float *cy, E_Float *cz)
{
  E_Int sizeNFace = cn.getSizeNFace();
  E_Int ncells = cn.getNElts();

  E_Int *count_neis = (E_Int *)XMALLOC(ncells * sizeof(E_Int));

  E_Float *lsqG = (E_Float *)XMALLOC(3*sizeNFace * sizeof(E_Float));
  E_Float *lsqGG = (E_Float *)XMALLOC(9*ncells * sizeof(E_Float));

  make_A_grad_matrices(cn, count_neis, owner, neigh, cx, cy, cz, lsqG);
  
  deduce_lsq_grad_matrices(cn, count_neis, lsqG, lsqGG);

  E_Float *b = (E_Float *)XMALLOC(sizeNFace * sizeof(E_Float));

  make_grad_RHS_vector(cn, count_neis, field, b, owner, neigh);

  E_Float *G = (E_Float *)XMALLOC(ncells*3 * sizeof(E_Float));
  make_gradient(cn, count_neis, b, G, lsqG, lsqGG); 

  XFREE(b);
  XFREE(count_neis);
  XFREE(lsqG);
  XFREE(lsqGG);

  return G;
}

static
void make_A_hess_matrices(K_FLD::FldArrayI &cn, E_Int *count_neis,
  E_Int *owner, E_Int *neigh, E_Float *cx, E_Float *cy, E_Float *cz,
  E_Float *lsqH)
{
  E_Int ncells = cn.getNElts();
  memset(count_neis, 0, ncells*sizeof(E_Int));

  E_Int nfaces = cn.getNFaces();
  E_Int *indPH = cn.getIndPH();

  // Internal faces
  for (E_Int i = 0; i < nfaces; i++) {
    E_Int nei = neigh[i];
    if (nei == -1) continue;

    E_Int own = owner[i];

    E_Float d[3] = {cx[nei]-cx[own],
                    cy[nei]-cy[own],
                    cz[nei]-cz[own]};

    E_Float *lo = &lsqH[6*(indPH[own] + count_neis[own]++)];
    E_Float *ln = &lsqH[6*(indPH[nei] + count_neis[nei]++)];

    E_Float dx = d[0];
    E_Float dy = d[1];
    E_Float dz = d[2];

    lo[0] = ln[0] = dx*dx;
    lo[1] = ln[1] = 2.0*dx*dy;
    lo[2] = ln[2] = 2.0*dx*dz;
    lo[3] = ln[3] = dy*dy;
    lo[4] = ln[4] = 2.0*dy*dz;
    lo[5] = ln[5] = dz*dz;
  }

  // TODO(Imad): boundary contributions
  // TODO(Imad): proc contributions
}

static
void deduce_lsq_hess_matrices(K_FLD::FldArrayI &cn, E_Int *count_neis,
  E_Float *lsqH, E_Float *lsqHH)
{
  E_Int ncells = cn.getNElts();
  E_Int *indPH = cn.getIndPH();

  for (E_Int i = 0; i < ncells; i++) {
    E_Float *pHH = &lsqHH[36*i];
    E_Float *pH = &lsqH[6*indPH[i]];
    E_Int idx = 0;
    for (E_Int j = 0; j < 6; j++) {
      for (E_Int k = 0; k < 6; k++) {
        pHH[idx] = 0.0;
        for (E_Int l = 0; l < count_neis[i]; l++)
          pHH[idx] += pH[l*6+j] * pH[l*6+k];
        idx++;
      }
    }
  }
}

static
void make_hess_RHS_vector(K_FLD::FldArrayI &cn, E_Int *count_neis,
  E_Float *field, E_Float *G, E_Float *b, E_Int *owner, E_Int *neigh,
  E_Float *cx, E_Float *cy, E_Float *cz)
{
  E_Int ncells = cn.getNElts();
  memset(count_neis, 0, ncells*sizeof(E_Int));

  E_Int nfaces = cn.getNFaces();
  E_Int *indPH = cn.getIndPH();

  for (E_Int i = 0; i < nfaces; i++) {
    E_Int nei = neigh[i];
    if (nei == -1) continue;

    E_Int own = owner[i];

    E_Float *go = &G[3*own];
    E_Float *gn = &G[3*nei];

    E_Float dfld = field[nei] - field[own];
    E_Float dx = cx[nei] - cx[own];
    E_Float dy = cy[nei] - cy[own];
    E_Float dz = cz[nei] - cz[own];

    b[indPH[own] + count_neis[own]++] =
      2.*( dfld - (go[0]*dx + go[1]*dy + go[2]*dz));
    b[indPH[nei] + count_neis[nei]++] =
      2.*(-dfld + (gn[0]*dx + gn[1]*dy + gn[2]*dz));
  }

  // TODO(Imad): boundary contributions
  // TODO(Imad): proc contributions
}

static
void make_hessian(K_FLD::FldArrayI &cn, E_Int *count_neis,
  E_Float *b, E_Float *H, E_Float *lsqH, E_Float *lsqHH)
{
  E_Int ncells = cn.getNElts();
  E_Int *indPH = cn.getIndPH();

  E_Float B[6];

  for (E_Int i = 0; i < ncells; i++) {
    // Make B vector = tA.b
    E_Float *pH = &lsqH[6*indPH[i]];
    E_Float *pb = &b[indPH[i]];
    for (E_Int j = 0; j < 6; j++) {
      B[j] = 0.0;
      for (E_Int k = 0; k < count_neis[i]; k++)
        B[j] += pb[k] * pH[k*6 + j];
    }

    // Solve
    E_Int converged = K_LINEAR::BiCGStab(&lsqHH[36*i], B, 6, &H[6*i]);

    if (!converged) {
      fprintf(stderr, "make_hessian: cell %d not converged.\n", i);
      exit(1);
    }
  }
}

PyObject *K_XCORE::computeHessian(PyObject *self, PyObject *args)
{
  PyObject *ARRAY, *FIELD, *CX, *CY, *CZ, *OWN, *NEI;
  if (!PYPARSETUPLE_(args, OOO_ OOO_ O_, &ARRAY, &FIELD, &CX, &CY, &CZ, &OWN,
      &NEI)) {
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
  //assert(ret == 1 && nfld == 1 && size == ncells);
  K_NUMPY::getFromNumpyArray(CY, cy, size, nfld, true);
  //assert(ret == 1 && nfld == 1 && size == ncells);
  K_NUMPY::getFromNumpyArray(CZ, cz, size, nfld, true);
  //assert(ret == 1 && nfld == 1 && size == ncells);

  // Parse owner and neigh
  E_Int *owner, *neigh;
  //E_Int nfaces = cn->getNFaces();
  K_NUMPY::getFromNumpyArray(OWN, owner, size, nfld, true);
  //assert(ret == 1 && nfld == 1 && size == nfaces);
  K_NUMPY::getFromNumpyArray(NEI, neigh, size, nfld, true);
  //assert(ret == 1 && nfld == 1 && size == nfaces);

  // Compute field gradient
  E_Float *G = compute_gradient(*cn, field, owner, neigh, cx, cy, cz);

  E_Int *count_neis = (E_Int *)XMALLOC(ncells * sizeof(E_Int));
  E_Int sizeNFace = cn->getSizeNFace();

  E_Float *lsqH = (E_Float *)XMALLOC(6*sizeNFace * sizeof(E_Float));

  E_Float *lsqHH = (E_Float *)XMALLOC(36*ncells * sizeof(E_Float));

  make_A_hess_matrices(*cn, count_neis, owner, neigh, cx, cy, cz, lsqH);

  deduce_lsq_hess_matrices(*cn, count_neis, lsqH, lsqHH);

  E_Float *b = (E_Float *)XMALLOC(sizeNFace * sizeof(E_Float));

  make_hess_RHS_vector(*cn, count_neis, field, G, b, owner, neigh, cx, cy, cz);

  npy_intp dims[2];
  dims[0] = (npy_intp)(ncells*6);
  dims[1] = 1;

  PyArrayObject *H = (PyArrayObject *)PyArray_SimpleNew(1, dims, NPY_DOUBLE);
  E_Float *pH = (E_Float *)PyArray_DATA(H);

  make_hessian(*cn, count_neis, b, pH, lsqH, lsqHH);

  XFREE(G);
  XFREE(owner);
  XFREE(neigh);
  XFREE(lsqH);
  XFREE(lsqHH);
  XFREE(count_neis);
  XFREE(b);

  return (PyObject *)H;
}
