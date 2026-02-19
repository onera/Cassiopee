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
# include "Array/Array.h"
#include "Connect/connect.h"

//=============================================================================
// Calcul du volume de toutes les cellules et des surfaces des interfaces
// CAS NON STRUCTURE
// IN: cn: nb de noeuds par elemt
// IN: eltType: list of Basic Element names
// IN: coordx, coordy, coordz: coordonnees x, y, z des pts de la grille
// OUT: snx, sny, snz: normales aux facettes %x, %y, %z
// OUT: surf: aires des facettes
// OUT: vol: volume des cellules
//=============================================================================
void K_METRIC::compMetricUnstruct(
  K_FLD::FldArrayI& cn, const char* eltType,
  const E_Float* coordx, const E_Float* coordy, const E_Float* coordz,
  E_Float* snx, E_Float* sny, E_Float* snz, E_Float* surf, E_Float* vol
)
{
  // Pre-compute element and facet offsets
  E_Int nc = cn.getNConnect();
  std::vector<char*> eltTypes;
  K_ARRAY::extractVars(eltType, eltTypes);

  // Number of facets per element
  std::vector<E_Int> nfpe;
  E_Int ierr = K_CONNECT::getNFPE(nfpe, eltType, false);
  if (ierr != 0) return;

  E_Int ntotFacets = 0;
  std::vector<E_Int> nepc(nc+1), nfpc(nc+1);
  nepc[0] = 0; nfpc[0] = 0;

  for (E_Int ic = 0; ic < nc; ic++)
  {
    K_FLD::FldArrayI& cm = *(cn.getConnect(ic));
    E_Int nelts = cm.getSize();
    nepc[ic+1] = nepc[ic] + nelts;
    nfpc[ic+1] = nfpc[ic] + nfpe[ic]*nelts;  // number of facets per connectivity
    ntotFacets += nfpe[ic]*nelts;
  }

  for (size_t ic = 0; ic < eltTypes.size(); ic++) delete [] eltTypes[ic];

  // Compute center of facets
  K_FLD::FldArrayF xint(ntotFacets), yint(ntotFacets), zint(ntotFacets);
  compUnstructCenterInt(cn, eltType, coordx, coordy, coordz,
                        xint.begin(), yint.begin(), zint.begin());

  // Compute facet normals and areas
  compSurfUnstruct(cn, eltType, coordx, coordy, coordz,
                   snx, sny, snz, surf);
  
  // Compute volume of elements
  compUnstructVol(cn, eltType, xint.begin(), yint.begin(), zint.begin(),
                  snx, sny, snz, vol);
}