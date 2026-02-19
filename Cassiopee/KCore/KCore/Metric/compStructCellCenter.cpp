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
# include "metric.h"


//=============================================================================
// Compute barycenter of cells.
// This version does not use minOfIntvoid
// IN: im, im, km: Number of mesh vertices along %i, %j, %k
// IN: xt, yt, zt: Vertex coordinates
// OUT: bary: Cell centers
//=============================================================================
void K_METRIC::compStructCellCenter(
  const E_Int im, const E_Int jm, const E_Int km,
  const E_Float* xt, const E_Float* yt, const E_Float* zt,
  E_Float* bary
)
{
  E_Int imjm = im * jm;
  E_Int im1 = im - 1;
  E_Int jm1 = jm - 1;
  E_Int km1 = km - 1;
  E_Int im1jm1 = im1 * jm1;

  E_Int inci = 1;
  E_Int incj = im;
  E_Int inck = imjm;
  E_Int incij = 1 + im;
  E_Int incik = 1 + imjm;
  E_Int incjk = im + imjm;
  E_Int incijk = 1 + im + imjm;

  #pragma omp parallel
  {
    E_Int i, j, k, rem;
    E_Int n1, n2, n3, n4, n5, n6, n7, n8;

    #pragma omp for
    for (E_Int lv = 0; lv < im1jm1 * km1; lv++)
    {
      k = lv / im1jm1;
      rem = lv % im1jm1;
      j = rem / im1;
      i = rem % im1;

      n1 = i + j * im + k * imjm;
      n2 = n1 + incj;
      n3 = n1 + incjk;
      n4 = n1 + inck;
      n5 = n1 + inci;
      n6 = n1 + incij;
      n7 = n1 + incijk;
      n8 = n1 + incik;

      bary[lv*3+0] = K_CONST::ONE_EIGHTH * (
        xt[n1] + xt[n2] + xt[n3] + xt[n4] +
        xt[n5] + xt[n6] + xt[n7] + xt[n8]);

      bary[lv*3+1] = K_CONST::ONE_EIGHTH * (
        yt[n1] + yt[n2] + yt[n3] + yt[n4] +
        yt[n5] + yt[n6] + yt[n7] + yt[n8]);

      bary[lv*3+2] = K_CONST::ONE_EIGHTH * (
        zt[n1] + zt[n2] + zt[n3] + zt[n4] +
        zt[n5] + zt[n6] + zt[n7] + zt[n8]);
    }
  }
}
