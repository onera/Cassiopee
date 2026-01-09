/*    
    Copyright 2013-2026 ONERA.

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
# include "post.h"

// ============================================================================
// Etant donnes n champs definis aux noeuds d une grille 3D, 
// calcul des champs aux centres des interfaces de la grille  
// fint1 = 0.25*(fa+fb+fc+fd)
// ============================================================================
void K_POST::compIntField(
  const E_Int ni, const E_Int nj, const E_Int nk,
  const E_Int nfld, const E_Float* f,
  E_Float* fint1
)
{
  E_Int ninj = ni * nj;
  E_Int ni1 = ni - 1;
  E_Int nj1 = nj - 1;
  E_Int nk1 = nk - 1;
  E_Int ni1nj1 = ni1 * nj1;
  E_Int ninj1 = ni * nj1;
  E_Int ni1nj = ni1 * nj;

  E_Int inti = ninj1 * nk1;
  E_Int intj = ni1nj * nk1;
  E_Int intij = inti + intj;

  // increments en i
  E_Int incij = ni;
  E_Int incik = ninj;
  E_Int incijk = ni + ninj;

  // increments en j
  E_Int incjk = ninj;
  E_Int incji = 1;
  E_Int incjik = 1 + ninj;

  // increments en k
  E_Int incki = 1;
  E_Int inckj = ni;
  E_Int inckij = 1 + ni;

  #pragma omp parallel
  {
    E_Int li, l1, l2, l3, l4;
    
    // Interfaces en i
    #pragma omp for collapse(3)
    for (E_Int k = 0; k < nk - 1; k++)
    for (E_Int j = 0; j < nj - 1; j++)
    for (E_Int i = 0; i <= ni1; i++)
    {
      li = i + j * ni + k * ninj1;
      l1 = i + j * ni + k * ninj;
      l2 = l1 + incij;
      l3 = l1 + incijk;
      l4 = l1 + incik;

      for (E_Int eq = 0; eq < nfld; eq++)
      {
        fint1[li * nfld + eq] = K_CONST::ONE_FOURTH *(
          f[l1 * nfld + eq] + f[l2 * nfld + eq] +
          f[l3 * nfld + eq] + f[l4 * nfld + eq]
        );
      }
    }
  
    // Interfaces en j
    #pragma omp for collapse(3)
    for (E_Int k = 0; k < nk - 1; k++)
    for (E_Int j = 0; j <= nj1; j++)
    for (E_Int i = 0; i < ni - 1; i++)
    {
      E_Int lj = i + j * ni1 + k * ni1nj + inti;
      E_Int l1 = i + j * ni + k * ninj;
      E_Int l2 = l1 + incjk;
      E_Int l3 = l1 + incjik;
      E_Int l4 = l1 + incji;

      for (E_Int eq = 0; eq < nfld; eq++)
      {
        fint1[lj * nfld + eq] = K_CONST::ONE_FOURTH *
          (f[l1 * nfld + eq] + f[l2 * nfld + eq] +
          f[l3 * nfld + eq] + f[l4 * nfld + eq]);
      }
    }

    // Interfaces en k
    #pragma omp for collapse(3)
    for (E_Int k = 0; k <= nk1; k++)
    for (E_Int j = 0; j < nj - 1; j++)
    for (E_Int i = 0; i < ni - 1; i++)
    {
      E_Int lk = i + j * ni1 + k * ni1nj1 + intij;
      E_Int l1 = i + j * ni + k * ninj;
      E_Int l2 = l1 + incki;
      E_Int l3 = l1 + inckij;
      E_Int l4 = l1 + inckj;

      for (E_Int eq = 0; eq < nfld; eq++)
      {
        fint1[lk * nfld + eq] = K_CONST::ONE_FOURTH *
          (f[l1 * nfld + eq] + f[l2 * nfld + eq] +
          f[l3 * nfld + eq] + f[l4 * nfld + eq]);
      }
    }
  }
}

// ============================================================================
// Compute the values of the vector (f1,f2,f3) at interfaces
// ============================================================================
void K_POST::compIntFieldV(
  const E_Int ni, const E_Int nj, const E_Int nk,
  const E_Float* f1, const E_Float* f2, const E_Float* f3,
  E_Float* fint1, E_Float* fint2, E_Float* fint3
)
{
  E_Int ninj = ni * nj;
  E_Int ni1 = ni - 1;
  E_Int nj1 = nj - 1;
  E_Int nk1 = nk - 1;
  E_Int ni1nj1 = ni1 * nj1;
  E_Int ninj1 = ni * nj1;
  E_Int ni1nj = ni1 * nj;

  E_Int inti = ninj1 * nk1;
  E_Int intj = ni1nj * nk1;
  E_Int intij = inti + intj;

  // increments en i
  E_Int incij = ni;
  E_Int incik = ninj;
  E_Int incijk = ni + ninj;

  // increments en j
  E_Int incjk = ninj;
  E_Int incji = 1;
  E_Int incjik = 1 + ninj;

  // increments en k
  E_Int incki = 1;
  E_Int inckj = ni;
  E_Int inckij = 1 + ni;

  #pragma omp parallel
  {
    E_Int li, lj, lk, l1, l2, l3, l4;

    // Interfaces en i
    #pragma omp for collapse(3)
    for (E_Int k = 0; k < nk - 1; k++)
    for (E_Int j = 0; j < nj - 1; j++)
    for (E_Int i = 0; i <= ni1; i++)
    {
      li = i + j * ni + k * ninj1;
      l1 = i + j * ni + k * ninj;
      l2 = l1 + incij;
      l3 = l1 + incijk;
      l4 = l1 + incik;

      fint1[li] = K_CONST::ONE_FOURTH * (f1[l1] + f1[l2] + f1[l3] + f1[l4]);
      fint2[li] = K_CONST::ONE_FOURTH * (f2[l1] + f2[l2] + f2[l3] + f2[l4]);
      fint3[li] = K_CONST::ONE_FOURTH * (f3[l1] + f3[l2] + f3[l3] + f3[l4]);
    }

    // Interfaces en j
    #pragma omp for collapse(3)
    for (E_Int k = 0; k < nk - 1; k++)
    for (E_Int j = 0; j <= nj1; j++)
    for (E_Int i = 0; i < ni - 1; i++)
    {
      lj = i + j * ni1 + k * ni1nj + inti;
      l1 = i + j * ni + k * ninj;
      l2 = l1 + incjk;
      l3 = l1 + incjik;
      l4 = l1 + incji;

      fint1[lj] = K_CONST::ONE_FOURTH * (f1[l1] + f1[l2] + f1[l3] + f1[l4]);
      fint2[lj] = K_CONST::ONE_FOURTH * (f2[l1] + f2[l2] + f2[l3] + f2[l4]);
      fint3[lj] = K_CONST::ONE_FOURTH * (f3[l1] + f3[l2] + f3[l3] + f3[l4]);
    }

    // Interfaces en k
    #pragma omp for collapse(3)
    for (E_Int k = 0; k <= nk1; k++)
    for (E_Int j = 0; j < nj - 1; j++)
    for (E_Int i = 0; i < ni - 1; i++)
    {
      lk = i + j * ni1 + k * ni1nj1 + intij;
      l1 = i + j * ni + k * ninj;
      l2 = l1 + incki;
      l3 = l1 + inckij;
      l4 = l1 + inckj;

      fint1[lk] = K_CONST::ONE_FOURTH * (f1[l1] + f1[l2] + f1[l3] + f1[l4]);
      fint2[lk] = K_CONST::ONE_FOURTH * (f2[l1] + f2[l2] + f2[l3] + f2[l4]);
      fint3[lk] = K_CONST::ONE_FOURTH * (f3[l1] + f3[l2] + f3[l3] + f3[l4]);
    }
  }
}
