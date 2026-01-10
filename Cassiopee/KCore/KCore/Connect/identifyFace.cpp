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

using namespace std;
using namespace K_FLD;

//=============================================================================
/* Identify a face in a basic elements connectivity
   IN: vecteurs d'indices des noeuds de la face (start 0)
   OUT: no de l'element et no local de la face
*/
//=============================================================================
E_Int K_CONNECT::identifyFace(E_Int* inds, E_Int n, 
                              vector< vector<E_Int> >& cVF)
{
  vector<E_Int> commonFaces;
  vector<E_Int> save;
  E_Int ind, indi;
  //printf("input ");
  //for (E_Int i = 0; i < n; i++) printf("%d ", inds[i]);
  //printf("\n");
  commonFaces = cVF[inds[0]];

  for (E_Int i = 1; i < n; i++)
  {
    // indice du noeud de la face (start 0)
    ind = inds[i];
    //printf("%d %d\n", i, ind);

    // faces candidates
    vector<E_Int>& cand = cVF[ind];
    
    // common part
    save = commonFaces; // copy
    commonFaces.clear();
    for (size_t i1 = 0; i1 < save.size(); i1++)
    {
      indi = save[i1];
      for (size_t i2 = 0; i2 < cand.size(); i2++)
        if (indi == cand[i2]) commonFaces.push_back(indi);
    }
    
    if (commonFaces.size() == 0) return -1; // not found
  }
  
  if (commonFaces.size() != 1) return -1; // error: not a face
  else return commonFaces[0];
}
