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

#include "Connect/connect.h"
#include <vector>
#include <algorithm>

using namespace K_FLD;
using namespace std;

//=============================================================================
/* Change a FE (nf, 2) connectivity of the type:
   for each interface, store the left and right elements into
   a NFace connectivity storing for each element, the number of faces and
   the face indices.
 */
//=============================================================================
void K_CONNECT::connectFE2NFace(FldArrayI& cFE, FldArrayI& cNFace, E_Int& nelts)
{
  E_Int nf = cFE.getSize();
  E_Int* facesp1 = cFE.begin(1);
  E_Int* facesp2 = cFE.begin(2);

  // Trouve le nbre d'elements (le max)
  E_Int e1, e2;
  nelts = 0;
  for (E_Int i = 0; i < nf; i++)
  {
    e1 = facesp1[i]; e2 = facesp2[i];
    nelts = K_FUNC::E_max(nelts, e1, e2);
  }
  vector< vector<E_Int> > a(nelts); // pour chaque element stocke les faces

  for (E_Int i = 0; i < nf; i++)
  {
    e1 = facesp1[i]-1; e2 = facesp2[i]-1;
    //printf("%d %d\n", e1, e2);
    if (e1 >= 0) a[e1].push_back(i+1);
    if (e2 >= 0) a[e2].push_back(i+1);
  }

  // unique it
  for (E_Int i = 0; i < nelts; i++)
  {
    sort(a[i].begin(), a[i].end());
    a[i].erase(unique(a[i].begin(), a[i].end()), a[i].end());
    //printf("%d %d\n", i, a[i].size());
  }

  // compactage
  E_Int size = 0;
  for (E_Int i = 0; i < nelts; i++)
  {
    size += (a[i].size()+1);
  }

  cNFace.malloc(size);
  E_Int* cn = cNFace.begin();

  size = 0;
  for (E_Int i = 0; i < nelts; i++)
  {
    nf = a[i].size();
    cn[size] = nf;
    for (E_Int j = 1; j <= nf; j++)
    { cn[size+j] = a[i][j-1]; }
    size += nf+1;
  }
}
