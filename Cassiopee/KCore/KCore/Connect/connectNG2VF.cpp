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
# include <vector>
# include <algorithm>
using namespace K_FLD;
using namespace std;

//=============================================================================
/* Change a NGON connectivity cNG where are stored two connectivities:
  cFV=face/vertex connectivity 
  cEF=elt/face connectivity 
  IN: cNG: NGON connectivity
  OUT: cVF: vertex/face connectivity, doit etre dimensionne au nb de noeuds*/
//=============================================================================
void K_CONNECT::connectNG2VF(FldArrayI& cNG, vector< vector<E_Int> >& cVF)
{
  // Acces non universel sur le ptrs
  E_Int* ngon = cNG.getNGon();
  E_Int* indPG = cNG.getIndPG();
  // Acces universel nbre de faces
  E_Int nfaces = cNG.getNFaces();
  
  E_Int nvert, node;
  for (E_Int nof = 0; nof < nfaces; nof++)
  {
    // Acces universel face nof
    E_Int* face = cNG.getFace(nof, nvert, ngon, indPG);
    for (E_Int nov = 0; nov < nvert; nov++)
    {
      node = face[nov]-1;
      cVF[node].push_back(nof+1);
    }
  }
  
  vector<E_Int>::iterator it;  
  for (unsigned int nov = 0; nov < cVF.size(); nov++)
  {
    sort(cVF[nov].begin(), cVF[nov].end()); 
    it = unique(cVF[nov].begin(), cVF[nov].end());
    cVF[nov].resize( it - cVF[nov].begin() );
  }
}