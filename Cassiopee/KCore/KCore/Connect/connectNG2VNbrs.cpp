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
# include "stdio.h"
using namespace K_FLD;
using namespace std;

//=============================================================================
/*
  Change a NGON connectivity cNG where are stored two connectivities,
  to a Vertex-Vertex neighbours connectivity
  IN: cNG: NGON connectivity
  OUT: cVN: Vertex-Vertex neighbours connectivity. For each vertex, 
  give the indices of its neighbours.
  Les index des vertex commence a 1.
  Attention: cVN doit deja etre alloue au nombre de noeuds.
*/
//=============================================================================
void K_CONNECT::connectNG2VNbrs(FldArrayI& cNG, vector< vector<E_Int> >& cVN)
{
  // Acces non universel sur le ptrs
  E_Int* ngon = cNG.getNGon();
  E_Int* nface = cNG.getNFace();
  E_Int* indPG = cNG.getIndPG();
  E_Int* indPH = cNG.getIndPH();
  // Acces universel nbre d'elements
  E_Int nelts = cNG.getNElts();
  // nombre total de noeuds
  E_Int nv = cVN.size();

  // construction de la connectivite Face/Noeuds (les noeuds sont definis plusieurs fois)
  E_Int nf, nv2, vertex1, vertex2;
  for (E_Int et = 0; et < nelts; et++)
  {
    // Acces universel element et
    E_Int* elt = cNG.getElt(et, nf, nface, indPH);
    for (E_Int j = 0; j < nf; j++)
    {
      // Acces universel face elt[j]-1
      E_Int* face = cNG.getFace(elt[j]-1, nv2, ngon, indPG);
      for (E_Int nov = 0; nov < nv2-1; nov++)
      {
        vertex1 = face[nov];
        vertex2 = face[nov+1];
        cVN[vertex1-1].push_back(vertex2);
        cVN[vertex2-1].push_back(vertex1);
      }
      // pour cycler
      vertex1 = face[nv2-1];
      vertex2 = face[0];
      cVN[vertex1-1].push_back(vertex2);
      cVN[vertex2-1].push_back(vertex1);
    }    
  } 

  // classement et unicite des noeuds voisins
  for (E_Int i = 0; i < nv; i++)
  {
    sort(cVN[i].begin(), cVN[i].end());
    cVN[i].erase(unique(cVN[i].begin(), cVN[i].end()), cVN[i].end());
  }
}