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
void compute_lsq_grad_matrices
(
  K_FLD::FldArrayI &cn, E_Float *x, E_Float *y, E_Float *z,
  E_Int *owner, E_Int *neigh, E_Int *count_neis, E_Float *centers,
  E_Float *lsqG, E_Float *lsqGG
)
{
  E_Int ncells = cn.getNElts();
  memset(count_neis, 0, ncells*sizeof(E_Int));

  // Internal faces
  E_Int nfaces = cn.getNFaces();
  E_Int *indPH = cn.getIndPH();
  assert(indPH[ncells] == 6*ncells);

  for (E_Int i = 0; i < nfaces; i++) {
    E_Int nei = neigh[i];
    if (nei == -1) continue;

    E_Int own = owner[i];
    
    E_Float *co = &centers[3*own];
    E_Float *cn = &centers[3*nei];

    E_Float d[3];
    for (E_Int j = 0; j < 3; j++)
      d[j] = cn[j] - co[j];

    E_Float *lo = &lsqG[3*(indPH[own] + count_neis[own]++)];
    E_Float *ln = &lsqG[3*(indPH[nei] + count_neis[nei]++)];

    for (E_Int j = 0; j < 3; j++) {
      lo[j] = d[j];
      ln[j] = -d[j];
    }
  }

  // TODO(Imad): boundary faces contribution

  // Compute tA.A
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
E_Int make_gradient
(
  K_FLD::FldArrayI &cn, E_Int *owner, E_Int *neigh, E_Int *count_neis,
  E_Float *centers, E_Float *fld, E_Float *b, E_Float *B, E_Float *lsqG,
  E_Float *lsqGG, E_Float *G
)
{
  E_Int ncells = cn.getNElts();
  memset(count_neis, 0, ncells*sizeof(E_Int));

  // Construct b vector
  E_Int nfaces = cn.getNFaces();
  E_Int *indPH = cn.getIndPH();
 
  for (E_Int i = 0; i < nfaces; i++) {
    E_Int nei = neigh[i];
    if (nei == -1) continue;

    E_Int own = owner[i];
    b[indPH[own] + count_neis[own]++] = fld[nei] - fld[own];
    b[indPH[nei] + count_neis[nei]++] = fld[own] - fld[nei];
  }

  // TODO(Imad): boundary faces contributions
  
  E_Int bad_gradient = 0;

  for (E_Int i = 0; i < ncells; i++) {
    // Construct B vector
    E_Float *pG = &lsqG[3*indPH[i]];
    E_Float *pb = &b[indPH[i]];
    for (E_Int j = 0; j < 3; j++) {
      B[j] = 0.0;
      for (E_Int k = 0; k < count_neis[i]; k++)
        B[j] += pb[k] * pG[k*3+j];
    }

    // Solve
    E_Int converged = BiCGStab(&lsqGG[9*i], &G[3*i], B, 3);
    bad_gradient |= !converged;
    if (!converged) {
      printf("Gradient: BiCGStab not converged\n");
      printf("Cell: " SF_D_ "\n", i);
      E_Float *pt = &lsqGG[9*i];
      for (E_Int j = 0; j < 3; j++) {
        for (E_Int k = 0; k < 3; k++) {
          printf("%.3e ", pt[3*j+k]);
        }
        puts("");
      }
      for (E_Int j = 0; j < 3; j++) {
        printf("%.3e ", B[j]);
      }
      puts("");
    }
  }

  return bad_gradient;
}

E_Int compute_gradients_ngon
(
  K_FLD::FldArrayI &cn, E_Float *x, E_Float *y, E_Float *z,
  E_Int *owner, E_Int *neigh, E_Float *centers,
  const std::vector<E_Float *> &flds,
  std::vector<E_Float *> &Gs
)
{
  E_Int ncells = cn.getNElts();
  E_Int *indPH = cn.getIndPH();

  E_Float *lsqG = (E_Float *)XMALLOC(3*indPH[ncells] * sizeof(E_Float));
  E_Float *lsqGG = (E_Float *)XMALLOC(9*ncells * sizeof(E_Float));
  E_Int *count_neis = (E_Int *)XMALLOC(ncells * sizeof(E_Int));
  
  compute_lsq_grad_matrices(cn, x, y, z, owner, neigh, count_neis,
    centers, lsqG, lsqGG);

  E_Float *B = (E_Float *)XMALLOC(3 * sizeof(E_Float));
  E_Float *b = (E_Float *)XMALLOC(indPH[ncells] * sizeof(E_Float));

  E_Int bad_gradient = 0;

  for (size_t i = 0; i < flds.size(); i++) {
    const auto &fld = flds[i];
    auto &G = Gs[i];
    bad_gradient |= make_gradient(cn, owner, neigh, count_neis, centers,
      fld, b, B, lsqG, lsqGG, G);
  }

  XFREE(lsqG);
  XFREE(lsqGG);
  XFREE(B);
  XFREE(b);
  XFREE(count_neis);

  return bad_gradient;
}
