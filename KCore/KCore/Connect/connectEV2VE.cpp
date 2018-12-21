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

using namespace K_FLD;
using namespace std;

//=============================================================================
/*
  Change a Elts-Vertex connectivity to a Vertex-Elts connectivity
  IN: cEV: Elts-Vertex connectivity. For each elt, give vertices index.
  OUT: cVE: Vertex-Elts connectivity. For each vertex, give elts index.
  Les indices de vertex commencent a 1
  Les indices d'elements commencent a 0
  Attention: cVE doit deja etre alloue au nombre de noeuds.
*/
//=============================================================================
void K_CONNECT::connectEV2VE(FldArrayI& cEV,
                             vector< vector<E_Int> >& cVE)
{
  E_Int ne = cEV.getSize();
  E_Int ns = cEV.getNfld();
  E_Int nv = cVE.size();
  E_Int vertex;

  // Calcul du nombre d'elements attaches a chaque noeud
  FldArrayI nve(nv); nve.setAllValuesAtNull();
  E_Int* nvep = nve.begin();

//#pragma omp parallel default (shared)
  {
    for (E_Int j = 1; j <= ns; j++)
    {
      E_Int* cEVp = cEV.begin(j);
//#pragma omp for
      for (E_Int i = 0; i < ne; i++)
      {
        vertex = cEVp[i]-1;
//#pragma omp atomic update
        nvep[vertex]++;
      }
    }
  }

#pragma omp parallel for default (shared)
  for (E_Int i = 0; i < nv; i++) cVE[i].reserve(nvep[i]);

  for (E_Int j = 1; j <= ns; j++)
  {
    E_Int* cEVp = cEV.begin(j);
    for (E_Int i = 0; i < ne; i++)
    {
      vertex = cEVp[i]-1;
      cVE[vertex].push_back(i);
    }
  }
}
