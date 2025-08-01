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
# include "post.h"
# include "Array/Array.h"

/* Calcul du gradient d'un champ defini aux noeuds d une grille non structuree
   retourne le gradient defini aux centres des elts.
   CAS 3D.
   IN: xt, yt, zt: coordonnees x, y, z des pts de la grille
   IN: cn: connectivite elts-noeuds
   IN: eltType: list of BE element types forming the ME mesh
   IN: field: champ defini aux noeuds auquel on applique grad
   OUT: fieldf: ieme champ defini aux facettes des elts
   OUT: snx, sny, snz: normales aux facettes %x, %y, %z
   OUT: surf: aires des facettes
   OUT: vol: Volume of the elements
   OUT: xint, yint, zint: Coordonnees du centre des facettes
   OUT: gradx, grady, gradz: gradient de field %x, %y, %z
*/
void K_POST::compUnstrGrad(
  const E_Float* xt, const E_Float* yt, const E_Float* zt,
  K_FLD::FldArrayI& cn, const char* eltType, const E_Float* field,
  E_Float* fieldf,
  E_Float* snx, E_Float* sny, E_Float* snz, E_Float* surf, E_Float* vol,
  E_Float* xint, E_Float* yint, E_Float* zint,
  E_Float* gradx, E_Float* grady, E_Float* gradz
)
{
  // Compute surface and volume of elements
  K_METRIC::compUnstrMetric(
    xt, yt, zt, cn,
    xint, yint, zint, snx, sny, snz, surf, vol
  );

  // Compute field on face centers
  compUnstrNodes2Faces(cn, field, fieldf);

  // Compute gradient at element centers
  E_Int nc = cn.getNConnect();
  std::vector<char*> eltTypes;
  K_ARRAY::extractVars(eltType, eltTypes);
  
  // Pre-compute element and face offsets
  vector<E_Int> nepc(nc+1), nfpc(nc+1);
  nepc[0] = 0; nfpc[0] = 0;
  for (E_Int ic = 0; ic < nc; ic++)
  {
    K_FLD::FldArrayI& cm = *(cn.getConnect(ic));
    E_Int nelts = cm.getSize();
    E_Int nvpe = cm.getNfld();

    nepc[ic+1] = nepc[ic] + nelts;
    nfpc[ic+1] = nfpc[ic] + nvpe*nelts;
  }
  
  #pragma omp parallel
  {
    E_Int pos, elOffset, fcOffset;
    E_Float sumx, sumy, sumz, invvol;
    
    for (E_Int ic = 0; ic < nc; ic++)
    {
      K_FLD::FldArrayI& cm = *(cn.getConnect(ic));
      E_Int nelts = cm.getSize();
      E_Int nvpe = cm.getNfld();
      elOffset = nepc[ic];
      fcOffset = nfpc[ic];

      #pragma omp for
      for (E_Int i = 0; i < nelts; i++)
      {
        sumx = K_CONST::E_ZERO_FLOAT;
        sumy = K_CONST::E_ZERO_FLOAT;
        sumz = K_CONST::E_ZERO_FLOAT;
        invvol = K_CONST::ONE / K_FUNC::E_max(vol[i], K_CONST::E_MIN_VOL);

        for (E_Int fi = 0; fi < nvpe; fi++)
        {
          pos = fcOffset + i * nvpe + fi;
          sumx += snx[pos] * fieldf[pos];
          sumy += sny[pos] * fieldf[pos];
          sumz += snz[pos] * fieldf[pos];
        }

        pos = elOffset + i;
        gradx[pos] = sumx * invvol;
        grady[pos] = sumy * invvol;
        gradz[pos] = sumz * invvol;
      }
    }
  }
}