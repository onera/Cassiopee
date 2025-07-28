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
// Compute barycenter of tetra cells.
// IN: cn: Element-Node connectivity
// IN: xt, yt, zt: Vertex coordinates
// OUT: bary: Cell centers
//=============================================================================
void K_METRIC::compTetraCellCenter(
  K_FLD::FldArrayI& cn,
  const E_Float* xt, const E_Float* yt, const E_Float* zt,
  E_Float* bary
)
{
  K_FLD::FldArrayI& cm = *(cn.getConnect(0));
  E_Int nelts = cm.getSize();
  
  #pragma omp parallel
  {
    E_Int et, ind1, ind2, ind3, ind4;

    #pragma omp for
    for (et = 0; et < nelts; et++)
    {
      ind1 = cm(et, 1) - 1;
      ind2 = cm(et, 2) - 1;
      ind3 = cm(et, 3) - 1;
      ind4 = cm(et, 4) - 1;

      bary[3*et+0] = K_CONST::ONE_FOURTH * (
        xt[ind1] + xt[ind2] + xt[ind3] + xt[ind4]);

      bary[3*et+1] = K_CONST::ONE_FOURTH * (
        yt[ind1] + yt[ind2] + yt[ind3] + yt[ind4]);

      bary[3*et+2] = K_CONST::ONE_FOURTH * (
        zt[ind1] + zt[ind2] + zt[ind3] + zt[ind4]);
    }
  }
}
