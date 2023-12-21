#include "Proto.h"

static
void make_gradient(AMesh *M, E_Int *count_neis, E_Float *b, E_Float *G)
{
  E_Float B[3];

  for (E_Int i = 0; i < M->ncells; i++) {
    // Make B vector = tA.b
    E_Float *pG = &M->lsqG[3*M->indPH[i]];
    E_Float *pb = &b[M->indPH[i]];
    for (E_Int j = 0; j < 3; j++) {
      B[j] = 0.0;
      for (E_Int k = 0; k < count_neis[i]; k++)
        B[j] += pb[k] * pG[k*3 + j];
    }

    // Solve
    E_Int converged = K_LINEAR::BiCGStab(&M->lsqGG[9*i], B, 3, &G[3*i]);

    assert(converged);
  }
}

static
void make_A_grad_matrices(AMesh *M, E_Int *count_neis)
{
  memset(count_neis, 0, M->ncells*sizeof(E_Int));

  // Internal faces
  for (E_Int i = 0; i < M->nfaces; i++) {
    E_Int nei = M->neigh[i];
    if (nei == -1) continue;

    E_Int own = M->owner[i];
    
    E_Float d[3] = {M->cx[nei]-M->cx[own],
                    M->cy[nei]-M->cy[own],
                    M->cz[nei]-M->cz[own]};

    E_Float *lo = &M->lsqG[3*(M->indPH[own] + count_neis[own]++)];
    E_Float *ln = &M->lsqG[3*(M->indPH[nei] + count_neis[nei]++)];

    for (E_Int j = 0; j < 3; j++) {
      lo[j] = d[j];
      ln[j] = -d[j];
    }
  }

  // TODO(Imad): boundary contributions

  // Proc faces
  for (E_Int i = 0; i < M->npatches; i++) {
    E_Int nf = M->patches[i].nfaces;
    E_Int *faces = M->patches[i].faces;
    E_Float *xn = M->xn[i];
    E_Float *yn = M->yn[i];
    E_Float *zn = M->zn[i];

    for (E_Int j = 0; j < nf; j++) {
      E_Int face = faces[j];
      E_Int own = M->owner[face];
      
      E_Float *lo = &M->lsqG[3*(M->indPH[own] + count_neis[own]++)];
      
      lo[0] = xn[j] - M->cx[own];
      lo[1] = yn[j] - M->cy[own];
      lo[2] = zn[j] - M->cz[own];
    }
  }
}

static
void deduce_lsq_grad_matrices(AMesh *M, E_Int *count_neis)
{
  for (E_Int i = 0; i < M->ncells; i++) {
    E_Float *pGG = &M->lsqGG[9*i];
    E_Float *pG = &M->lsqG[3*M->indPH[i]];
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
void make_grad_RHS_vector(AMesh *M, E_Int *count_neis, E_Float *field,
  E_Float *b)
{
  memset(count_neis, 0, M->ncells*sizeof(E_Int));

  // Construct b vector
  for (E_Int i = 0; i < M->nfaces; i++) {
    E_Int nei = M->neigh[i];
    if (nei == -1) continue;

    E_Int own = M->owner[i];
    b[M->indPH[own] + count_neis[own]++] = field[nei] - field[own];
    b[M->indPH[nei] + count_neis[nei]++] = field[own] - field[nei];
  }

  // TODO(Imad): boundary faces contributions

  // Proc faces
  for (E_Int i = 0; i < M->npatches; i++) {
    E_Int nf = M->patches[i].nfaces;
    E_Int *faces = M->patches[i].faces;
    E_Float *fn = M->fn[i];

    for (E_Int j = 0; j < nf; j++) {
      E_Int face = faces[j];
      E_Int own = M->owner[face];
      b[M->indPH[own] + count_neis[own]++] = fn[j] - field[own];
    }
  }
}


static
E_Float *compute_gradient(AMesh *M, E_Float *field)
{
  E_Int sizeNFace = M->indPH[M->ncells];

  E_Int *count_neis = (E_Int *)XMALLOC(M->ncells * sizeof(E_Int));
 
  // Compute lsq matrices only once per adaptation run
  if (!M->lsqG) {
    M->lsqG = (E_Float *)XMALLOC(3*sizeNFace * sizeof(E_Float));
    M->lsqGG = (E_Float *)XMALLOC(9*M->ncells * sizeof(E_Float));

    make_A_grad_matrices(M, count_neis);
  
    deduce_lsq_grad_matrices(M, count_neis);
  }

  E_Float *b = (E_Float *)XMALLOC(sizeNFace * sizeof(E_Float));

  make_grad_RHS_vector(M, count_neis, field, b);

  E_Float *G = (E_Float *)XMALLOC(M->ncells*3 * sizeof(E_Float));
  make_gradient(M, count_neis, b, G); 

  XFREE(b);
  XFREE(count_neis);

  return G;
}

static
void make_A_hess_matrices(AMesh *M, E_Int *count_neis)
{
  memset(count_neis, 0, M->ncells*sizeof(E_Int));

  // Internal faces
  for (E_Int i = 0; i < M->nfaces; i++) {
    E_Int nei = M->neigh[i];
    if (nei == -1) continue;

    E_Int own = M->owner[i];

    E_Float d[3] = {M->cx[nei]-M->cx[own],
                    M->cy[nei]-M->cy[own],
                    M->cz[nei]-M->cz[own]};

    E_Float *lo = &M->lsqH[6*(M->indPH[own] + count_neis[own]++)];
    E_Float *ln = &M->lsqH[6*(M->indPH[nei] + count_neis[nei]++)];

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

  // Proc faces
  for (E_Int i = 0; i < M->npatches; i++) {
    E_Int nf = M->patches[i].nfaces;
    E_Int *faces = M->patches[i].faces;
    E_Float *xn = M->xn[i];
    E_Float *yn = M->yn[i];
    E_Float *zn = M->zn[i];

    for (E_Int j = 0; j < nf; j++) {
      E_Int own = M->owner[faces[j]];

      E_Float d[3] = {xn[j] - M->cx[own],
                      yn[j] - M->cy[own],
                      zn[j] - M->cz[own]};

      E_Float *lo = &M->lsqH[6*(M->indPH[own] + count_neis[own]++)];

      E_Float dx = d[0];
      E_Float dy = d[1];
      E_Float dz = d[2];

      lo[0] = dx*dx;
      lo[1] = 2.0*dx*dy;
      lo[2] = 2.0*dx*dz;
      lo[3] = dy*dy;
      lo[4] = 2.0*dy*dz;
      lo[5] = dz*dz;
    }
  }
}

static
void deduce_lsq_hess_matrices(AMesh *M, E_Int *count_neis)
{
  for (E_Int i = 0; i < M->ncells; i++) {
    E_Float *pHH = &M->lsqHH[36*i];
    E_Float *pH = &M->lsqH[6*M->indPH[i]];
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
void make_hess_RHS_vector(AMesh *M, E_Int *count_neis, E_Float *field,
  E_Float *G, E_Float *b)
{
  memset(count_neis, 0, M->ncells*sizeof(E_Int));

  for (E_Int i = 0; i < M->nfaces; i++) {
    E_Int nei = M->neigh[i];
    if (nei == -1) continue;

    E_Int own = M->owner[i];

    E_Float *go = &G[3*own];
    E_Float *gn = &G[3*nei];

    E_Float dfld = field[nei] - field[own];
    E_Float dx = M->cx[nei] - M->cx[own];
    E_Float dy = M->cy[nei] - M->cy[own];
    E_Float dz = M->cz[nei] - M->cz[own];

    b[M->indPH[own] + count_neis[own]++] =
      2.*( dfld - (go[0]*dx + go[1]*dy + go[2]*dz));
    b[M->indPH[nei] + count_neis[nei]++] =
      2.*(-dfld + (gn[0]*dx + gn[1]*dy + gn[2]*dz));
  }

  // TODO(Imad): boundary contributions

  // Proc faces
  for (E_Int i = 0; i < M->npatches; i++) {
    E_Int nf = M->patches[i].nfaces;
    E_Int *faces = M->patches[i].faces;
    E_Float *fld = M->fn[i];
    E_Float *xn = M->xn[i];
    E_Float *yn = M->yn[i];
    E_Float *zn = M->zn[i];

    for (E_Int j = 0; j < nf; j++) {
      E_Int own = M->owner[faces[j]];

      E_Float *go = &G[3*own];

      E_Float dfld = fld[j] - field[own];
      E_Float dx = xn[j] - M->cx[own];
      E_Float dy = yn[j] - M->cy[own];
      E_Float dz = zn[j] - M->cz[own];

      b[M->indPH[own] + count_neis[own]++] =
        2.*(dfld - (go[0]*dx + go[1]*dy + go[2]*dz));
    }
  }
}

static
void make_hessian(AMesh *M, E_Int *count_neis, E_Float *b, E_Float *H)
{
  E_Float B[6];

  for (E_Int i = 0; i < M->ncells; i++) {
    // Make B vector = tA.b
    E_Float *pH = &M->lsqH[6*M->indPH[i]];
    E_Float *pb = &b[M->indPH[i]];
    for (E_Int j = 0; j < 6; j++) {
      B[j] = 0.0;
      for (E_Int k = 0; k < count_neis[i]; k++)
        B[j] += pb[k] * pH[k*6 + j];
    }

    // Solve
    E_Int converged = K_LINEAR::BiCGStab(&M->lsqHH[36*i], B, 6, &H[6*i]);

    assert(converged);
  }
}

E_Float *compute_hessian(AMesh *M, E_Float *field)
{
  exchange_proc_data_d(M, field, &M->fn);

  E_Float *G = compute_gradient(M, field);

  E_Int *count_neis = (E_Int *)XMALLOC(M->ncells * sizeof(E_Int));

  E_Int sizeNFace = M->indPH[M->ncells];

  if (!M->lsqH) {
    M->lsqH = (E_Float *)XMALLOC(6*sizeNFace * sizeof(E_Float *));
    M->lsqHH = (E_Float *)XMALLOC(36*M->ncells * sizeof(E_Float *));
    
    make_A_hess_matrices(M, count_neis);

    deduce_lsq_hess_matrices(M, count_neis);
  }

  E_Float *b = (E_Float *)XMALLOC(sizeNFace * sizeof(E_Float));

  make_hess_RHS_vector(M, count_neis, field, G, b);

  E_Float *H = (E_Float *)XMALLOC(6*M->ncells * sizeof(E_Float));
  make_hessian(M, count_neis, b, H); 

  XFREE(G);
  XFREE(count_neis);
  XFREE(b);

  return H;
}