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
// Compute cell barycenter for a ME mesh.
// IN: cn: Element-Node connectivity
// IN: xt, yt, zt: Vertex coordinates
// OUT: xb, yb, zb: Barycenter coordinates
//=============================================================================
void K_METRIC::compUnstructCellCenter(
  K_FLD::FldArrayI& cn,
  const E_Float* xt, const E_Float* yt, const E_Float* zt,
  E_Float* xb, E_Float* yb, E_Float* zb
)
{
  // Pre-compute the number of element per connectivity (gives the offsets)
  E_Int nc = cn.getNConnect();
  std::vector<E_Int> nepc(nc+1);
  nepc[0] = 0;
  for (E_Int ic = 0; ic < nc; ic++)
  {
    K_FLD::FldArrayI& cm = *(cn.getConnect(ic));
    E_Int nelts = cm.getSize();
    nepc[ic+1] = nepc[ic] + nelts;
  }

  #pragma omp parallel
  {
    E_Int ind, pos, nelts, nvpe;

    for (E_Int ic = 0; ic < nc; ic++)
    {
      K_FLD::FldArrayI& cm = *(cn.getConnect(ic));
      nelts = cm.getSize();
      nvpe = cm.getNfld();
  
      #pragma omp for
      for (E_Int i = 0; i < nelts; i++)
      {
        pos = nepc[ic] + i;
        for (E_Int j = 1; j <= nvpe; j++)
        {
          ind = cm(i, j) - 1;
          xb[pos] += xt[ind] / nvpe;
          yb[pos] += yt[ind] / nvpe;
          zb[pos] += zt[ind] / nvpe;
        }
      }
    }
  }
}
