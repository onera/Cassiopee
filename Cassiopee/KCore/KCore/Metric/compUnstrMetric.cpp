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

//=============================================================================
// Calcul du volume de toutes les cellules et des surfaces des interfaces
// CAS NON STRUCTURE
// IN: cn: nb de noeuds par elemt
// IN: eltType: list of Basic Element names
// IN: coordx, coordy, coordz: coordonnees x, y, z des pts de la grille
// OUT: xint, yint, zint: coordonnees du centre des facettes 
// OUT: snx, sny, snz: normales aux facettes %x, %y, %z
// OUT: surf: aires des facettes
// OUT: vol: volume des cellules
//=============================================================================
void K_METRIC::compUnstructMetric(
  K_FLD::FldArrayI& cn, const char* eltType,
  const E_Float* coordx, const E_Float* coordy, const E_Float* coordz,
  E_Float* xint, E_Float* yint, E_Float* zint,
  E_Float* snx, E_Float* sny, E_Float* snz,
  E_Float* surf, E_Float* vol
)
{
  // Compute center of facets
  compUnstructCenterInt(cn, eltType, coordx, coordy, coordz,
                        xint, yint, zint);

  // Compute facet normals and areas
  compUnstructSurf(cn, eltType, coordx, coordy, coordz,
                   snx, sny, snz, surf);

  // Compute volume of elements
  E_Int fcOffset = 0;
  E_Int elOffset = 0;
  E_Int nc = cn.getNConnect();
  std::vector<char*> eltTypes;
  K_ARRAY::extractVars(eltType, eltTypes);

  for (E_Int ic = 0; ic < nc; ic++)
  {
    K_FLD::FldArrayI& cm = *(cn.getConnect(ic));
    E_Int nelts = cm.getSize();
    E_Int nfpe = -1;            // number of facets per element
    E_Int nfpc = nfpe*nelts;    // number of facets per connectivity

    if (strcmp(eltTypes[ic], "TRI") == 0) nfpe = 1;
    else if (strcmp(eltTypes[ic], "QUAD") == 0) nfpe = 1;
    else if (strcmp(eltTypes[ic], "TETRA") == 0) nfpe = 4;
    else if (strcmp(eltTypes[ic], "HEXA") == 0) nfpe = 6;
    else if (strcmp(eltTypes[ic], "PENTA") == 0) nfpe = 5;
    else if (strcmp(eltTypes[ic], "PYRA") == 0) nfpe = 5;
    else
    {
      fprintf(stderr, "Error: in K_METRIC::compUnstructMetric.\n");
      fprintf(stderr, "Unknown type of element.\n");
      exit(0);
    }

    #pragma omp parallel
    {
      E_Int pos;
      E_Float voli;

      #pragma omp for
      for (E_Int i = 0; i < nelts; i++)
      {
        voli = K_CONST::E_ZERO_FLOAT;
        for (E_Int fidx = 0; fidx < nfpe; fidx++)
        {
          pos = fcOffset + i * nfpe + fidx;
          voli += xint[pos] * snx[pos]
                + yint[pos] * sny[pos]
                + zint[pos] * snz[pos];
        }
        vol[elOffset+i] = K_CONST::ONE_THIRD * voli;
      }
    }

    // Increment offsets
    fcOffset += nfpc;
    elOffset += nelts;
  }
}