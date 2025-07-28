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
// Calcul les vecteurs normaux aux triangles
// IN: cn: connectivite elts-noeuds
// IN: xt, yt, zt: coordonnees x, y, z des pts de la grille
// OUT: nsurf: normales aux facettes
//=============================================================================
void K_METRIC::compNormUnstructSurf(
  K_FLD::FldArrayI& cn,
  const E_Float* xt, const E_Float* yt, const E_Float* zt,
  E_Float* nsurf
)
{
  K_FLD::FldArrayI& cm = *(cn.getConnect(0));
  E_Int nt = cm.getSize();

  #pragma omp parallel
  {
    E_Int ind1, ind2, ind3;
    E_Float l1x, l1y, l1z;
    E_Float l2x, l2y, l2z;
    E_Float surfx, surfy, surfz;

    #pragma omp for
    for (E_Int i = 0; i < nt; i++)
    {
      ind1 = cm(i, 1);
      ind2 = cm(i, 2);
      ind3 = cm(i, 3);

      l1x = xt[ind1] - xt[ind2];
      l1y = yt[ind1] - yt[ind2];
      l1z = zt[ind1] - zt[ind2];

      l2x = xt[ind1] - xt[ind3];
      l2y = yt[ind1] - yt[ind3];
      l2z = zt[ind1] - zt[ind3];

      surfx = l1y*l2z - l1z*l2y;
      surfy = l1z*l2x - l1x*l2z;
      surfz = l1x*l2y - l1y*l2x;

      nsurf[i*3+0] = K_CONST::ONE_HALF * surfx;
      nsurf[i*3+1] = K_CONST::ONE_HALF * surfy;
      nsurf[i*3+2] = K_CONST::ONE_HALF * surfz;
    }
  }
}

void K_METRIC::compNormUnstructSurft(
  K_FLD::FldArrayI& cn,
  const E_Float* xt, const E_Float* yt, const E_Float* zt,
  E_Float* nxt, E_Float* nyt, E_Float* nzt
)
{
  K_FLD::FldArrayI& cm = *(cn.getConnect(0));
  E_Int nt = cm.getSize();
  
  #pragma omp parallel
  {
    E_Int ind1, ind2, ind3;
    E_Float l1x, l1y, l1z;
    E_Float l2x, l2y, l2z;
    E_Float surfx, surfy, surfz;

    #pragma omp for
    for (E_Int i = 0; i < nt; i++)
    {
      ind1 = cm(i, 1);
      ind2 = cm(i, 2);
      ind3 = cm(i, 3);

      l1x = xt[ind1] - xt[ind2];
      l1y = yt[ind1] - yt[ind2];
      l1z = zt[ind1] - zt[ind2];

      l2x = xt[ind1] - xt[ind3];
      l2y = yt[ind1] - yt[ind3];
      l2z = zt[ind1] - zt[ind3];

      surfx = l1y * l2z - l1z * l2y;
      surfy = l1z * l2x - l1x * l2z;
      surfz = l1x * l2y - l1y * l2x;

      nxt[i] = K_CONST::ONE_HALF * surfx;
      nyt[i] = K_CONST::ONE_HALF * surfy;
      nzt[i] = K_CONST::ONE_HALF * surfz;
    }
  }
}
