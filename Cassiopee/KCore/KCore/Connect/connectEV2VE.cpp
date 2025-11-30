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
void K_CONNECT::connectEV2VE(FldArrayI& cEV, vector< vector<E_Int> >& cVE)
{
  // Acces universel sur BE/ME
  E_Int nc = cEV.getNConnect();

  E_Int npts = cVE.size(); // Nombre de points du maillage
  // Calcul du nombre d'elements attaches a chaque noeud
  FldArrayI nve(npts); nve.setAllValuesAtNull();
  E_Int* nvep = nve.begin();
  E_Int vertex;
  E_Int offset = 0; // decalage

  // Boucle sur toutes les connectivites pour pre-calculer la valence de
  // chaque noeud
  for (E_Int ic = 0; ic < nc; ic++)
  {
    FldArrayI& cm = *(cEV.getConnect(ic));
    // Nbre elements de cette connectivite
    E_Int nelts = cm.getSize();
    // Nombre de points par elements de cette connectivite
    E_Int nvpe = cm.getNfld();

    for (E_Int j = 1; j <= nvpe; j++)
    {
      for (E_Int i = 0; i < nelts; i++)
      {
        vertex = cm(i, j);
        nvep[vertex-1]++;
      }
    }
  }

#pragma omp parallel for default (shared)
  for (E_Int i = 0; i < npts; i++) cVE[i].reserve(nvep[i]);

  // Boucle sur toutes les connectivites pour remplir cVE
  for (E_Int ic = 0; ic < nc; ic++)
  {
    FldArrayI& cm = *(cEV.getConnect(ic));
    E_Int nelts = cm.getSize();
    E_Int nvpe = cm.getNfld();
    
    for (E_Int j = 1; j <= nvpe; j++)
    {
      for (E_Int i = 0; i < nelts; i++)
      {
        vertex = cm(i, j);
        cVE[vertex-1].push_back(offset + i);
      }
    }

    offset += nelts; // increment element offset
  }
}
