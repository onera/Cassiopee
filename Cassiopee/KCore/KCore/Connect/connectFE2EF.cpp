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

using namespace K_FLD;

//=============================================================================
/* 
   Pour les NGONS conformes.
   IN: cFE: connectivite faces->elements (nfaces,2)
   IN: nelts: nbre d'elements
   OUT: cEF: connectivite elements->faces (comme dans la connectivite NG
   complete: Npts de la face, indices des pts, ..., Npts de la face, ...
   Cette routine dimensionne cEF.
*/
//=============================================================================
void K_CONNECT::connectFE2EF(FldArrayI& cFE, E_Int nelts, FldArrayI& cEF)
{
  E_Int nfaces = cFE.getSize();
  // Compte les faces pour chaque elements
  FldArrayI nbface(nelts); nbface.setAllValuesAtNull();
  E_Int* nbfacep = nbface.begin();
  E_Int a, b;
  E_Int* cn1 = cFE.begin(1);
  E_Int* cn2 = cFE.begin(2);
  for (E_Int i = 0; i < nfaces; i++)
  {
    a = cn1[i]; b = cn2[i];
    if (a != 0) nbfacep[a-1]++;
    if (b != 0) nbfacep[b-1]++;
  }
  
  // position des faces pour chaque elt
  FldArrayI pos(nelts);
  E_Int* posp = pos.begin();
  posp[0] = 0;
  for (E_Int i = 0; i < nelts-1; i++) posp[i+1] = posp[i] + nbfacep[i] + 1;
  
  E_Int sizeEF = nelts;
  for (E_Int i = 0; i < nelts; i++) sizeEF += nbfacep[i];
  cEF.malloc(sizeEF); 
  E_Int* cEFp = cEF.begin();
  
  #pragma omp parallel
  {
    #pragma omp for
    for (E_Int i = 0; i < nelts; i++) cEFp[posp[i]] = nbfacep[i];
    #pragma omp for
    for (E_Int i = 0; i < nelts; i++) nbfacep[i] = 0;
  }
  
  E_Int ind1, ind2;
  for (E_Int i = 0; i < nfaces; i++)
  {
    a = cn1[i]-1; b = cn2[i]-1;
    if (a+1 != 0)
    {
      ind1 = posp[a] + nbfacep[a] + 1;
      cEFp[ind1] = i+1; nbfacep[a]++;
    }
    if (b+1 != 0)
    {
      ind2 = posp[b] + nbfacep[b] + 1;
      cEFp[ind2] = i+1; nbfacep[b]++;
    }
  }
}
