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

#define MAXNEI 6 // ok for hexa, tetra and prism

static
void compute_lsq_hess_matrices
(
  K_FLD::FldArrayI &cn, E_Float *x, E_Float *y, E_Float *z,
  E_Int *owner, E_Int *neigh, E_Int *count_neis, E_Float *centers,
  E_Float *lsqH, E_Float *lsqHH
)
{
  E_Int ncells = cn.getNElts();
  memset(count_neis, 0, ncells*sizeof(E_Int));

  // Internal faces
  E_Int nfaces = cn.getNFaces();
  E_Int *indPH = cn.getIndPH();

  for (E_Int i = 0; i < nfaces; i++) {
    E_Int nei = neigh[i];
    if (nei == -1) continue;

    E_Int own = owner[i];
    
    E_Float *co = &centers[3*own];
    E_Float *cn = &centers[3*nei];
    
    E_Float dx = cn[0] - co[0];
    E_Float dy = cn[1] - co[1];
    E_Float dz = cn[2] - co[2];
    
    E_Float *lo = &lsqH[6*(indPH[own] + count_neis[own]++)];
    E_Float *ln = &lsqH[6*(indPH[nei] + count_neis[nei]++)];

    lo[0] = dx*dx;
    lo[1] = 2.*dx*dy;
    lo[2] = 2.*dx*dz;
    lo[3] = dy*dy;
    lo[4] = 2.*dy*dz;
    lo[5] = dz*dz;

    ln[0] = dx*dx;
    ln[1] = 2.*dx*dy;
    ln[2] = 2.*dx*dz;
    ln[3] = dy*dy;
    ln[4] = 2.*dy*dz;
    ln[5] = dz*dz;
  }

  // TODO(Imad): boundary faces contribution

  // Compute tA.A
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
E_Int make_hessian
(
  K_FLD::FldArrayI &cn, E_Int *owner, E_Int *neigh, E_Int *count_neis,
  E_Float *centers, E_Float *G, E_Float *fld, E_Float *b, E_Float *B,
  E_Float *lsqH, E_Float *lsqHH, E_Float *H
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

    E_Float *co = &centers[3*own];
    E_Float *cn = &centers[3*nei];

    E_Float *go = &G[3*own];
    E_Float *gn = &G[3*nei];

    E_Float dfld = fld[nei] - fld[own];
    E_Float dx = cn[0] - co[0];
    E_Float dy = cn[1] - co[1];
    E_Float dz = cn[2] - co[2];

    b[indPH[own] + count_neis[own]++] =
      2.*( dfld - (go[0]*dx + go[1]*dy + go[2]*dz));
    b[indPH[nei] + count_neis[nei]++] =
      2.*(-dfld + (gn[0]*dx + gn[1]*dy + gn[2]*dz));
  }

  // TODO(Imad): boundary faces contributions
  E_Int bad_hessian = 0;
  for (E_Int i = 0; i < ncells; i++) {
    // construct B vector
    E_Float *pH = &lsqH[6*indPH[i]];
    E_Float *pb = &b[indPH[i]];
    for (E_Int j = 0; j < 6; j++) {
      B[j] = 0.0;
      for (E_Int k = 0; k < count_neis[i]; k++)
        B[j] += pb[k] * pH[k*6+j];
    }

    // solve
    E_Int converged = BiCGStab(&lsqHH[36*i], &H[6*i], B, 6);
    bad_hessian |= !converged;
    if (!converged) {
      printf("BiCGStab not converged\n");
      printf("Cell: " SF_D_ "\n", i);
      E_Float *pt = &lsqHH[36*i];
      for (E_Int j = 0; j < 6; j++) {
        for (E_Int k = 0; k < 6; k++) {
          printf("%.3e ", pt[6*j+k]);
        }
        puts("");
      }
      for (E_Int j = 0; j < 6; j++) {
        printf("%.3e ", B[j]);
      }
      puts("");
    }
  }

  return bad_hessian;
}

E_Int compute_hessians_ngon
(
  K_FLD::FldArrayI &cn, E_Float *x, E_Float *y, E_Float *z,
  E_Int *owner, E_Int *neigh, E_Float *centers,
  const std::vector<E_Float *> &Gs,
  const std::vector<E_Float *> &flds,
  std::vector<E_Float *> &Hs
)
{
  E_Int ncells = cn.getNElts();
  E_Int *indPH = cn.getIndPH();

  E_Float *lsqH = (E_Float *)XMALLOC(6*indPH[ncells] * sizeof(E_Float));
  E_Float *lsqHH = (E_Float *)XMALLOC(36*ncells * sizeof(E_Float));
  E_Int *count_neis = (E_Int *)XMALLOC(ncells * sizeof(E_Int));

  compute_lsq_hess_matrices(cn, x, y, z, owner, neigh, count_neis, 
    centers, lsqH, lsqHH);

  E_Float *B = (E_Float *)XMALLOC(6 * sizeof(E_Float));
  E_Float *b = (E_Float *)XMALLOC(indPH[ncells] * sizeof(E_Float));

  E_Int bad_hessian = 0;

  for (size_t i = 0; i < flds.size(); i++) {
    const auto &fld = flds[i];
    const auto &G = Gs[i];
    auto &H = Hs[i];
    bad_hessian |= make_hessian(cn, owner, neigh, count_neis, centers, G,
      fld, b, B, lsqH, lsqHH, H);
  }

  XFREE(lsqH);
  XFREE(lsqHH);
  XFREE(B);
  XFREE(b);
  XFREE(count_neis);

  return bad_hessian;
}
