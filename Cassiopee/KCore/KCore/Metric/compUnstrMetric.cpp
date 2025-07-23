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
// Calcul du volume de toutes les cellules et des surfaces des interfaces
// CAS NON STRUCTURE
// IN: npts: nb de pts du maillage
// IN: nelts: nb d elements
// IN: nedges: nb de facettes par elemt
// IN: nnodes: nb de noeuds par elemt
// IN: cn: nb de noeuds par elemt
// IN: coordx, coordy, coordz: coordonnees x, y, z des pts de la grille
// OUT: xint, yint, zint: coordonnees du centre des facettes 
// OUT: snx, sny, snz: normales aux facettes %x, %y, %z
// OUT: surf: aires des facettes
// OUT: vol: volume des cellules
//=============================================================================
void K_METRIC::compUnstructMetric(
  const E_Int nedges, const E_Int nnodes,
  K_FLD::FldArrayI& cn,
  const E_Float* coordx, const E_Float* coordy, const E_Float* coordz,
  E_Float* xint, E_Float* yint, E_Float* zint,
  E_Float* snx, E_Float* sny, E_Float* snz,
  E_Float* surf, E_Float* vol
)
{
  K_FLD::FldArrayI& cm = *(cn.getConnect(0));
  E_Int nelts = cm.getSize();
  
  // Compute center of facets
  compUnstructCenterInt(nedges, cn,
                        coordx, coordy, coordz, xint, yint, zint);

  // Compute facet normals and surfaces
  compUnstructSurf(nedges, cn,
                   coordx, coordy, coordz, snx, sny, snz, surf);

  #pragma omp parallel
  {
    E_Float voli;
    #pragma omp for
    for (E_Int i = 0; i < nelts; i++)
    {
      voli = K_CONST::E_ZERO_FLOAT;
      for (E_Int edge = 0; edge < nedges; edge++)
      {
        voli += xint[i*nedges+edge] * snx[i*nedges+edge]
              + yint[i*nedges+edge] * sny[i*nedges+edge]
              + zint[i*nedges+edge] * snz[i*nedges+edge];
      }
      vol[i] = K_CONST::ONE_THIRD * voli;
    }
  }
}