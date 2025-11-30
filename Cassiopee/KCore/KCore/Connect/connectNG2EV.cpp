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
# include "Connect/connect.h"

//=============================================================================
// Fill an Element-Vertex connectivity from an NGON connectivity.
//   IN: cNG: NGON connectivity
//   OUT: cEV: Element-Vertex connectivity: vertex numbering starts at 1, vertex
//             indices are sorted
//   cEV must have previously been allocated
//=============================================================================
void K_CONNECT::connectNG2EV(
  K_FLD::FldArrayI& cNG,
  std::vector<std::vector<E_Int> >& cEV
)
{
  E_Int* ngon = cNG.getNGon(); E_Int* nface = cNG.getNFace();
  E_Int* indPG = cNG.getIndPG(); E_Int* indPH = cNG.getIndPH();
  E_Int nelts = cNG.getNElts();

  #pragma omp parallel
  {
    E_Int nf, nv, fidx;
    
    #pragma omp for
    for (E_Int i = 0; i < nelts; i++)
    {
      // Loop over each face of this element and store all of their vertices
      std::vector<E_Int>& vertices = cEV[i];
      vertices.clear(); vertices.reserve(24);
      E_Int* elt = cNG.getElt(i, nf, nface, indPH);
      for (E_Int j = 0; j < nf; j++)
      {
        fidx = elt[j] - 1;
        E_Int* face = cNG.getFace(fidx, nv, ngon, indPG);
        for (E_Int v = 0; v < nv; v++) vertices.push_back(face[v]);
      }

      // Remove duplicated vertices
      std::sort(vertices.begin(), vertices.end());
      vertices.erase(
        std::unique(vertices.begin(), vertices.end()), vertices.end()
      );
    }
  }
}