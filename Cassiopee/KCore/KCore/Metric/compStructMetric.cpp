/*
    Copyright 2013-2025 Onera.

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
// Calcul du volume de toutes les cellules et des surface des interfaces.
// CAS STRUCTURE
// IN: ni, nj, nk: Number of mesh vertices along %i, %j, %k
// IN: nbInti, nbIntj, nbIntk: Number of interfaces along the %i, %j, %k directions
// IN: x, y, z: Vertex coordinates
// OUT: vol: Volume of the elements
// OUT: surfx, surfy, surfz: Surface vectors
// OUT: snorm: Norm of the surface vectors
// OUT: cix, ciy, ciz: Interface centers
//=============================================================================
void K_METRIC::compStructMetric(
  const E_Int ni, const E_Int nj, const E_Int nk,
  const E_Int nbInti, const E_Int nbIntj, const E_Int nbIntk,
  const E_Float* x, const E_Float* y, const E_Float* z,
  E_Float* vol, E_Float* surfx, E_Float* surfy, E_Float* surfz,
  E_Float* snorm, E_Float* cix, E_Float* ciy, E_Float* ciz
)
{
  E_Int nbIntij, ni1, nj1, nk1, ni1nj1;
  E_Int incInti, incIntj, incIntk;

  ni1 = ni - 1;
  nj1 = nj - 1;
  nk1 = nk - 1;
  ni1nj1 = ni1 * nj1;
  nbIntij = nbInti + nbIntj;

  compCenterInterface(ni, nj, nk, nbInti, nbIntij, x, y, z,
                      cix, ciy, ciz);

  compIntSurf(ni, nj, nk, nbInti, nbIntij, x, y, z,
              surfx, surfy, surfz, snorm);

  incInti = 1;
  incIntj = ni1;
  incIntk = ni1nj1;

  #pragma omp parallel
  {
    E_Int i, j, k, lv;
    E_Int n1, n2, n3, n4, n5, n6;
    E_Float v1, v2, v3, v4, v5, v6;
    
    #pragma omp for
    for (E_Int c = 0; c < ni1nj1 * nk1; c++)
    {
      i = c % ni1;
      j = (c / ni1) % nj1;
      k = c / ni1nj1;

      lv = i + j * ni1 + k * ni1nj1;
      n1 = i + j * ni + k * ni * nj1;
      v1 = cix[n1] * surfx[n1] + ciy[n1] * surfy[n1] + ciz[n1] * surfz[n1];

      n2 = n1 + incInti;
      v2 = cix[n2] * surfx[n2] + ciy[n2] * surfy[n2] + ciz[n2] * surfz[n2];

      n3 = i + j * ni1 + k * ni1 * nj + nbInti;
      v3 = cix[n3] * surfx[n3] + ciy[n3] * surfy[n3] + ciz[n3] * surfz[n3];

      n4 = n3 + incIntj;
      v4 = cix[n4] * surfx[n4] + ciy[n4] * surfy[n4] + ciz[n4] * surfz[n4];

      n5 = i + j * ni1 + k * ni1nj1 + nbIntij;
      v5 = cix[n5] * surfx[n5] + ciy[n5] * surfy[n5] + ciz[n5] * surfz[n5];

      n6 = n5 + incIntk;
      v6 = cix[n6] * surfx[n6] + ciy[n6] * surfy[n6] + ciz[n6] * surfz[n6];

      vol[lv] = K_CONST::ONE_THIRD * (v2 - v1 + v4 - v3 + v6 - v5);
    }
  }
}
