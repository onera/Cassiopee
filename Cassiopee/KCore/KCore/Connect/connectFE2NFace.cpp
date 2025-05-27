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

// NGON(CGNSv3)
void K_CONNECT::connectFE2NFace(FldArrayI& cFE, FldArrayI& cNFace, E_Int& nelts)
{
  E_Int nf = cFE.getSize();
  E_Int* facesp1 = cFE.begin(1);
  E_Int* facesp2 = cFE.begin(2);

  // Trouve le nbre d'elements (le max)
  nelts = 0;
  #pragma omp parallel for reduction(max:nelts)
  for (E_Int i = 0; i < nf; i++)
  {
    E_Int e1 = facesp1[i]; E_Int e2 = facesp2[i];
    nelts = K_FUNC::E_max(nelts, e1, e2);
  }
  vector< vector<E_Int> > a(nelts); // stocke les faces pour chaque element
  for (E_Int i = 0; i < nelts; i++) a[i].reserve(6);

  for (E_Int i = 0; i < nf; i++)
  {
    E_Int e1 = facesp1[i]-1; E_Int e2 = facesp2[i]-1;
    if (e1 >= 0) a[e1].push_back(i+1);
    if (e2 >= 0) a[e2].push_back(i+1);
  }

  // unique it
  for (E_Int i = 0; i < nelts; i++)
  {
    sort(a[i].begin(), a[i].end());
    a[i].erase(unique(a[i].begin(), a[i].end()), a[i].end());
  }

  // compactage
  E_Int size = 0;
  FldArrayI off(nelts);
  E_Int* o = off.begin();
  for (E_Int i = 0; i < nelts; i++)
  {
    o[i] = size;
    size += a[i].size()+1;
  }

  cNFace.malloc(size);
  E_Int* cn = cNFace.begin();

  #pragma omp parallel
  {
    E_Int oi;
    #pragma omp for
    for (E_Int i = 0; i < nelts; i++)
    {
      nf = a[i].size(); oi = o[i];
      cn[oi] = nf;
      for (E_Int j = 1; j <= nf; j++) cn[oi+j] = a[i][j-1];
    }
  }
  
  off.malloc(0);
}

// NGON(CGNSv3) + index
void K_CONNECT::connectFE2NFace3(FldArrayI& cFE, FldArrayI& cNFace, FldArrayI& off, E_Int& nelts)
{
  E_Int nf = cFE.getSize();
  E_Int* facesp1 = cFE.begin(1);
  E_Int* facesp2 = cFE.begin(2);

  // Trouve le nbre d'elements (le max)
  nelts = 0;
  #pragma omp parallel for reduction(max:nelts)
  for (E_Int i = 0; i < nf; i++)
  {
    E_Int e1 = facesp1[i]; E_Int e2 = facesp2[i];
    nelts = K_FUNC::E_max(nelts, e1, e2);
  }
  vector< vector<E_Int> > a(nelts); // stocke les faces pour chaque element
  for (E_Int i = 0; i < nelts; i++) a[i].reserve(6);

  for (E_Int i = 0; i < nf; i++)
  {
    E_Int e1 = facesp1[i]-1; E_Int e2 = facesp2[i]-1;
    if (e1 >= 0) a[e1].push_back(i+1);
    if (e2 >= 0) a[e2].push_back(i+1);
  }

  // unique it
  for (E_Int i = 0; i < nelts; i++)
  {
    sort(a[i].begin(), a[i].end());
    a[i].erase(unique(a[i].begin(), a[i].end()), a[i].end());
  }

  // compactage
  E_Int size = 0;
  off.malloc(nelts);
  E_Int* o = off.begin();
  for (E_Int i = 0; i < nelts; i++)
  {
    o[i] = size;
    size += a[i].size()+1;
  }

  cNFace.malloc(size);
  E_Int* cn = cNFace.begin();

  #pragma omp parallel
  {
    E_Int oi;
    #pragma omp for
    for (E_Int i = 0; i < nelts; i++)
    {
      nf = a[i].size(); oi = o[i];
      cn[oi] = nf;
      for (E_Int j = 1; j <= nf; j++) cn[oi+j] = a[i][j-1];
    }
  }
}

// NGON(CGNSv4) + offset
void K_CONNECT::connectFE2NFace4(FldArrayI& cFE, FldArrayI& cNFace, FldArrayI& off, E_Int& nelts)
{
  E_Int nf = cFE.getSize();
  E_Int* facesp1 = cFE.begin(1);
  E_Int* facesp2 = cFE.begin(2);

  // Trouve le nbre d'elements (le max)
  nelts = 0;
  #pragma omp parallel for reduction(max:nelts)
  for (E_Int i = 0; i < nf; i++)
  {
    E_Int e1 = facesp1[i]; E_Int e2 = facesp2[i];
    nelts = K_FUNC::E_max(nelts, e1, e2);
  }
  vector< vector<E_Int> > a(nelts); // stocke les faces pour chaque element
  for (E_Int i = 0; i < nelts; i++) a[i].reserve(6);

  for (E_Int i = 0; i < nf; i++)
  {
    E_Int e1 = facesp1[i]-1; E_Int e2 = facesp2[i]-1;
    if (e1 >= 0) a[e1].push_back(i+1);
    if (e2 >= 0) a[e2].push_back(i+1);
  }

  // unique it
  for (E_Int i = 0; i < nelts; i++)
  {
    sort(a[i].begin(), a[i].end());
    a[i].erase(unique(a[i].begin(), a[i].end()), a[i].end());
  }

  // compactage
  E_Int size = 0;
  off.malloc(nelts+1);
  E_Int* o = off.begin();
  for (E_Int i = 0; i < nelts; i++)
  {
    o[i] = size;
    size += a[i].size();
  }
  o[nelts] = size;

  cNFace.malloc(size);
  E_Int* cn = cNFace.begin();

  #pragma omp parallel
  {
    E_Int oi;
    #pragma omp for
    for (E_Int i = 0; i < nelts; i++)
    {
      nf = a[i].size(); oi = o[i];
      for (E_Int j = 0; j < nf; j++) cn[oi+j] = a[i][j];
    }
  }
}
