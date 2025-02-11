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
# include <vector>
# include <algorithm>
using namespace K_FLD;
using namespace std;

//=============================================================================
/* Change a NGON connectivity cNG where are stored two connectivities:
  cFV=face/vertex connectivity 
  cEF=elt/face connectivity 
  IN: cNG: NGON connectivity
  OUT: cEV: elt/vertex connectivity: numerotation des vertices demarre a 1 
  cEV doit etre alloue au nb d'elements */
//=============================================================================
void K_CONNECT::connectNG2EV(FldArrayI& cNG, vector< vector<E_Int> >& cEV)
{
  // Acces non universel sur le ptrs
  E_Int* ngon = cNG.getNGon();
  E_Int* nface = cNG.getNFace();
  E_Int* indPG = cNG.getIndPG();
  E_Int* indPH = cNG.getIndPH();
  // Acces universel nbre d'elements
  E_Int ncells = cNG.getNElts();

  #pragma omp parallel
  {
    E_Int nf, nv;
    
    #pragma omp for
    for (E_Int et = 0; et < ncells; et++)
    {
      vector<E_Int>& vertices = cEV[et]; // noeuds associes a l'element et
      // Acces universel element et
      E_Int* elt = cNG.getElt(et, nf, nface, indPH);
      for (E_Int j = 0; j < nf; j++)
      {
        // Acces universel face elt[j]-1
        E_Int* face = cNG.getFace(elt[j]-1, nv, ngon, indPG);
        for (E_Int nov = 0; nov < nv; nov++)
        {
          vertices.push_back(face[nov]);
        }
      }
    }
  }
  vector<E_Int>::iterator it;
  for (E_Int et = 0; et < ncells; et++)
  {
    sort(cEV[et].begin(), cEV[et].end()); 
    it = unique(cEV[et].begin(), cEV[et].end()); 
    cEV[et].resize(it - cEV[et].begin());
  }
}