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
#include "Connect/connect.h"

//=============================================================================
/* Fill a Face-Element connectivity from an NGON connectivity.
  IN: cNG: NGON connectivity (conformal mesh)
  OUT: cFE: Face-Element connectivity. If a face has 1 neighbour element only,
            the second index value is 0.
  Size of cFE set in this routine.
 */
//=============================================================================
void K_CONNECT::connectNG2FE(K_FLD::FldArrayI& cNG, K_FLD::FldArrayI& cFE)
{
  E_Int* nface = cNG.getNFace(); E_Int* indPH = cNG.getIndPH();
  E_Int nfaces = cNG.getNFaces(); E_Int nelts = cNG.getNElts();

  // Face-Element connectivity
  cFE.malloc(nfaces, 2); cFE.setAllValuesAtNull();
  E_Int* facesp1 = cFE.begin(1);
  E_Int* facesp2 = cFE.begin(2);

  E_Int nf, fidx;
  for (E_Int i = 0; i < nelts; i++)
  {
    E_Int* elt = cNG.getElt(i, nf, nface, indPH);
    for (E_Int j = 0; j < nf; j++)
    {
      fidx = elt[j] - 1;
      if (facesp1[fidx] == 0) facesp1[fidx] = i+1;
      else facesp2[fidx] = i+1;
    }
  }
}
