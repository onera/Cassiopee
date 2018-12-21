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
/* Calcule la connectivite Elts/Elts voisins a partir de la connectivite 
   Face/Elts pour un NGON conforme.
   cEEN doit etre deja dimensionne au nombre d'elements. */
//=============================================================================
void K_CONNECT::connectFE2EENbrs(FldArrayI& cFE, 
                                 vector< vector<E_Int> >& cEEN)
{
  E_Int nf = cFE.getSize();
  
  E_Int e1, e2;
  E_Int* cFE1 = cFE.begin(1);
  E_Int* cFE2 = cFE.begin(2);
  for (E_Int i = 0; i < nf; i++)
  {
    e1 = cFE1[i]; e2 = cFE2[i];
    if (e1 > 0 && e2 > 0) {
      cEEN[e1-1].push_back(e2-1);
      cEEN[e2-1].push_back(e1-1); }
  }
}
