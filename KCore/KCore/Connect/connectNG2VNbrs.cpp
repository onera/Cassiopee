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
  E_Int* cnp = cNG.begin();             // pointeur sur la connectivite
  E_Int nfaces = cnp[0];
  E_Int sizeFN = cnp[1];                // dimension du tableau de connectivite Face/Noeuds
  E_Int nelts = cnp[sizeFN+2];          // nombre d'elements
  vector< vector<E_Int> > cEV(nelts);   // connectivite Element/Noeuds
  E_Int nv = cVN.size();                // nombre total de noeuds
  E_Int vertex1, vertex2;
  FldArrayI posFaces(nfaces); // tableau de position des faces dans la connectivite
  K_CONNECT::getPosFaces(cNG, posFaces);
  E_Int* posFacesp = posFaces.begin(); // pointeur sur posFace
  E_Int pos; // position d une face donnee dans la connectivite

  // 1- construction de la connectivite Face/Noeuds (les noeuds sont definis plusieurs fois)
  E_Int* ptr = cnp+sizeFN+4;//debut connectivite EF
  E_Int nf, face, nv2;
  for (E_Int et = 0; et < nelts; et++)
  {
    nf = ptr[0]; //nb de faces pour l elt
    for (E_Int j = 1; j <= nf; j++)
    {
      face = ptr[j]-1;// numero de la face
      pos = posFacesp[face];
      nv2 = cnp[pos];//nb de noeuds pour la face courante 
      for (E_Int nov = 1; nov < nv2; nov++)
      {
        vertex1 = cnp[pos+nov];//demarre a 1
        vertex2 = cnp[pos+nov+1];//demarre a 1
        cVN[vertex1-1].push_back(vertex2);
        cVN[vertex2-1].push_back(vertex1);
      }
      // pour cycler
      vertex1 = cnp[pos+nv2];//demarre a 1
      vertex2 = cnp[pos+1];//demarre a 1
      cVN[vertex1-1].push_back(vertex2);
      cVN[vertex2-1].push_back(vertex1);
    }    
    ptr += nf+1;
  }  

  // classement et unicite des noeuds voisins
  for (E_Int i = 0; i < nv; i++)
  {
    sort(cVN[i].begin(), cVN[i].end());
    cVN[i].erase(unique(cVN[i].begin(), cVN[i].end()), cVN[i].end());
  }
}
