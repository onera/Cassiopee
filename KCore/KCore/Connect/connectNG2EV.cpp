/*    
    Copyright 2013-2019 Onera.

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
  E_Int* cnp = cNG.begin();
  E_Int nfaces = cnp[0];
  E_Int sizeFN = cnp[1];
  E_Int ncells = cnp[sizeFN+2];
  E_Int* ptr = cnp+sizeFN+4;//debut connectivite EF
  E_Int nf, nv, vert;
  FldArrayI posFaces(nfaces); // tableau de position des faces dans la connectivite
  K_CONNECT::getPosFaces(cNG, posFaces);
  E_Int* posFacesp = posFaces.begin(); // pointeur sur posFace
  E_Int pos; // position d'une face donnee dans la connectivite

  for (E_Int et = 0; et < ncells; et++)
  {
    vector<E_Int>& vertices = cEV[et]; // noeuds associes a l'element et
    nf = ptr[0]; //nb de faces pour l elt
    for (E_Int j = 1; j <= nf; j++)
    {
      pos = posFacesp[ptr[j]-1];
      
      nv = cnp[pos];//nb de noeuds pour la face courante 
      for (E_Int nov = 1; nov <= nv; nov++)
      {
        vert = cnp[pos+nov];//demarre a 1
        vertices.push_back(vert);
      }
    }    
    ptr+= nf+1;
  }  
  vector<E_Int>::iterator it;  
  for (E_Int et = 0; et < ncells; et++)
  {
    sort(cEV[et].begin(), cEV[et].end()); 
    it = unique(cEV[et].begin(), cEV[et].end()); 
    cEV[et].resize(it - cEV[et].begin());
  }
}
