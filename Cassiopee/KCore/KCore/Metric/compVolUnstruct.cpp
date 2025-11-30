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
# include "Array/Array.h"
#include "Connect/connect.h"

//=============================================================================
// Calcul du volume des elements pour un maillage Multi-Elements.
// IN: cn: Element-Node connectivity
// IN: eltType: Element names
// IN: xint, yint, zint: coordonnees du centre des facettes
// IN: snx, sny, snz: normales aux facettes %x, %y, %z
// OUT: vol: volume des cellules
//=============================================================================
void K_METRIC::compUnstructVol(
  K_FLD::FldArrayI& cn, const char* eltType,
  const E_Float* xint, const E_Float* yint, const E_Float* zint,
  const E_Float* snx, const E_Float* sny, const E_Float* snz, E_Float* vol
)
{
  E_Int nc = cn.getNConnect();
  std::vector<char*> eltTypes;
  K_ARRAY::extractVars(eltType, eltTypes);

  // Number of facets per element
  std::vector<E_Int> nfpe;
  E_Int ierr = K_CONNECT::getNFPE(nfpe, eltType, false);
  if (ierr != 0) return;

  std::vector<E_Int> nepc(nc+1), nfpc(nc+1);
  nepc[0] = 0; nfpc[0] = 0;

  for (E_Int ic = 0; ic < nc; ic++)
  {
    K_FLD::FldArrayI& cm = *(cn.getConnect(ic));
    E_Int nelts = cm.getSize();
    nepc[ic+1] = nepc[ic] + nelts;
    nfpc[ic+1] = nfpc[ic] + nfpe[ic]*nelts;
  }

  for (size_t ic = 0; ic < eltTypes.size(); ic++) delete [] eltTypes[ic];

  #pragma omp parallel
  {
    E_Int pos, nelts, elOffset, fctOffset;
    E_Float voli;

    for (E_Int ic = 0; ic < nc; ic++)
    {
      K_FLD::FldArrayI& cm = *(cn.getConnect(ic));
      nelts = cm.getSize();
      elOffset = nepc[ic];
      fctOffset = nfpc[ic];

      #pragma omp for
      for (E_Int i = 0; i < nelts; i++)
      {
        voli = K_CONST::E_ZERO_FLOAT;
        for (E_Int fidx = 0; fidx < nfpe[ic]; fidx++)
        {
          pos = fctOffset + i * nfpe[ic] + fidx;
          voli += xint[pos] * snx[pos]
                + yint[pos] * sny[pos]
                + zint[pos] * snz[pos];
        }
        vol[elOffset + i] = K_CONST::ONE_THIRD * voli;
      }
    }
  }
}