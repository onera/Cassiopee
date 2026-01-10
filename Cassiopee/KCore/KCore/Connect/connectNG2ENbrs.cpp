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

# include "Connect/connect.h"
# include <vector>
# include <stdlib.h>
using namespace K_FLD;
using namespace std;

//=============================================================================
/*
  Calcule les elements contenants le sommet d'indice indV pour un NGON conforme. 
  IN: indV: indice du point considere
  IN: cVF: connectivite Vertex/Faces
  IN: cFE: connectivite Faces/Elts
  OUT: ENbrs: indices des elements contenant le sommet indV
*/
//=============================================================================
void K_CONNECT::connectNG2ENbrs(E_Int indV, vector< vector<E_Int> >& cVF, FldArrayI& cFE, 
                                vector<E_Int>& ENbrs)
{
  // cFE est de la forme (nfaces, 2)
  E_Int* facep1 = cFE.begin(1);  //element d'un cote de la face
  E_Int* facep2 = cFE.begin(2);  //element de l'autre cote de la face (0 si sur le bord)
  
  for (size_t i = 0; i < cVF[indV].size(); i++)
  { 
    E_Int face = cVF[indV][i]-1;
    if (facep1[face] != 0) ENbrs.push_back(facep1[face]-1);
    if (facep2[face] != 0) ENbrs.push_back(facep2[face]-1);
  }
  sort(ENbrs.begin(), ENbrs.end());
  ENbrs.erase(unique(ENbrs.begin(), ENbrs.end()), ENbrs.end());
}
