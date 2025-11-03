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

#include "Connect/connect.h"

using namespace K_FLD;

//=============================================================================
/* Change a NGON connectivity cNG where are stored two connectivities:
  cFV=face/vertex connectivity
  cEF=elt/face connectivity
  into a face/element connectivity.
  IN: cNG: connectivite NGON conforme
  OUT: cFE: face/elts connectivity
  if a face has not 2 neighbouring elements, its value is 0.
  cFE est dimensionne par cette routine.
 */
//=============================================================================
void K_CONNECT::connectNG2FE(FldArrayI& cNG, FldArrayI& cFE)
{
  E_Int* nface = cNG.getNFace(); E_Int* indPH = cNG.getIndPH();
  E_Int nfaces = cNG.getNFaces(); E_Int nelts = cNG.getNElts();

  // Connectivite faces->elts
  cFE.malloc(nfaces, 2); cFE.setAllValuesAtNull();
  E_Int* facesp1 = cFE.begin(1);
  E_Int* facesp2 = cFE.begin(2);

  E_Int nf, fidx;
  for (E_Int i = 0; i < nelts; i++)
  {
    E_Int* elt = cNG.getElt(i, nf, nface, indPH);
    for (E_Int j = 0; j < nf; j++)
    {
      fidx = K_FUNC::E_abs(elt[j]) - 1;
      if (facesp1[fidx] == 0) facesp1[fidx] = i+1;
      else facesp2[fidx] = i+1;
    }
  }
}
