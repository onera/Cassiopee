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

// Calcul de la divergence d'un champ défini aux noeuds d'une grille non structurée
// Retourne la divergence définie aux centres des éléments
// CAS 3D
// ============================================================================

#include "Def/DefFortranConst.h"

void compUnstrDiv(
  const E_Int dim, const E_Int npts, const E_Int nelts,
  const E_Int nedges,
  const E_Int nnodes,
  const E_Int* cn,
  const E_Float* xt,
  const E_Float* yt,
  const E_Float* zt,
  const E_Float* fieldX,
  const E_Float* fieldY,
  const E_Float* fieldZ,
  E_Float* fieldfx,
  E_Float* fieldfy,
  E_Float* fieldfz,
  E_Float* snx,
  E_Float* sny,
  E_Float* snz,
  E_Float* surf,
  E_Float* vol,
  E_Float* xint,
  E_Float* yint,
  E_Float* zint,
  E_Float* div
)
{
  // local vars
  E_Int eti, fi;
  E_Float gradxx, gradyy, gradzz;
  E_Float invvol;

  // calcul de la surface + volume des elts
  compUnstrMetric(
    npts, nelts, nedges, nnodes, cn,
    xt, yt, zt, xint, yint, zint, snx, sny, snz, surf, vol
  );

  // calcul des champs aux centres des facettes
  unstrNodes2Faces(dim, npts, nelts, nedges, nnodes, cn, fieldX, fieldfx);
  unstrNodes2Faces(dim, npts, nelts, nedges, nnodes, cn, fieldY, fieldfy);
  unstrNodes2Faces(dim, npts, nelts, nedges, nnodes, cn, fieldZ, fieldfz);

  // calcul de la divergence au centre des elts
  #pragma omp parallel private(eti, gradxx, gradyy, gradzz, invvol, fi)
  {
    #pragma omp for
    for (eti = 0; eti < nelts; eti++)
    {
      gradxx = K_CONST::ZERO;
      gradyy = K_CONST::ZERO;
      gradzz = K_CONST::ZERO;

      for (fi = 0; fi < nedges; fi++)
      {
        gradxx = gradxx + snx[eti*nedges + fi] * fieldfx[eti*nedges + fi];
        gradyy = gradyy + sny[eti*nedges + fi] * fieldfy[eti*nedges + fi];
        gradzz = gradzz + snz[eti*nedges + fi] * fieldfz[eti*nedges + fi];
      }

      invvol = K_CONST::ONE / fmax(vol[eti], E_MIN_VOL);
      div[eti] = (gradxx + gradyy + gradzz) * invvol;
    }
  }
}
