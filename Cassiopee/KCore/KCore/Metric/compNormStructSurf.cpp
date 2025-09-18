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
// Attention cette routine est 2D: pas de calcul des normales sur un hexaedre!
// IN: ni, nj: Number of mesh vertices along %i, %j
// IN: xt, yt, zt: Vertex coordinates
// OUT: nxt, nyt, nzt: Norm of the surface vector
//=============================================================================
void K_METRIC::compNormStructSurf(
  const E_Int ni, const E_Int nj,
  const E_Float* xt, const E_Float* yt, const E_Float* zt,
  E_Float* nxt, E_Float* nyt, E_Float* nzt
)
{
  #pragma omp parallel
  {
    E_Int indcell, ind1, ind2, ind3, ind4;
    E_Float l1x, l1y, l1z;
    E_Float l2x, l2y, l2z;
    E_Float x1, y1, z1;
    E_Float surf1x, surf1y, surf1z;
    E_Float surf2x, surf2y, surf2z;

    #pragma omp for
    for (E_Int j = 0; j < nj - 1; j++)
    {
      for (E_Int i = 0; i < ni - 1; i++)
      {
        ind1 = i + j * ni;
        ind2 = ind1 + 1;
        ind3 = ind2 + ni;
        ind4 = ind1 + ni;
        indcell = i + j * (ni - 1);
        
        // AB x AC
        x1 = xt[ind1];
        y1 = yt[ind1];
        z1 = zt[ind1];

        l1x = xt[ind2] - x1;
        l1y = yt[ind2] - y1;
        l1z = zt[ind2] - z1;

        l2x = xt[ind3] - x1;
        l2y = yt[ind3] - y1;
        l2z = zt[ind3] - z1;

        surf1x = l1y * l2z - l1z * l2y;
        surf1y = l1z * l2x - l1x * l2z;
        surf1z = l1x * l2y - l1y * l2x;

        // AC x AD
        l1x = xt[ind3] - x1;
        l1y = yt[ind3] - y1;
        l1z = zt[ind3] - z1;

        l2x = xt[ind4] - x1;
        l2y = yt[ind4] - y1;
        l2z = zt[ind4] - z1;

        surf2x = l1y * l2z - l1z * l2y;
        surf2y = l1z * l2x - l1x * l2z;
        surf2z = l1x * l2y - l1y * l2x;

        nxt[indcell] = K_CONST::ONE_HALF * (surf1x + surf2x);
        nyt[indcell] = K_CONST::ONE_HALF * (surf1y + surf2y);
        nzt[indcell] = K_CONST::ONE_HALF * (surf1z + surf2z);
      }   
    }
  }
}