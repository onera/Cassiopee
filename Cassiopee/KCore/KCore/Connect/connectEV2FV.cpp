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

#include "kcore.h"
#include "Connect/connect.h"

using namespace K_FLD;
using namespace std;

//=============================================================================
/* Change un connectivite Elts-Vertex (basic elements) en une connectivite
   Faces->Vertex. L'indice des faces est global, soit : nof + nelt*nfaces
   ou nfaces est le nbre de faces de l'elements, nelt le no de l'element
   (commencant a 0) et nof la numero local de la face 0,1,2,...
*/
//=============================================================================
void K_CONNECT::connectEV2FV(FldArrayI& cEV, const char* eltType,
                             FldArrayI& cFV)
{
  // Acces universel sur BE/ME
  E_Int nc = cEV.getNConnect();
  // Acces universel aux eltTypes
  vector<char*> eltTypes;
  K_ARRAY::extractVars(eltType, eltTypes);
  // Number of elements and faces per connectivity
  vector<E_Int> nelts(nc);
  vector<E_Int> nfaces(nc);
  vector<E_Int> nof(nc);
  // Accumulated number of faces over connectivities
  vector<E_Int> ntotfaces(nc+1); ntotfaces[0] = 0;
  vector<vector<E_Int> > face(nc);
    
  // Boucle sur toutes les connectivites pour remplir face et determiner
  // la taille de cFV
  for (E_Int ic = 0; ic < nc; ic++)
  {
    FldArrayI& cm = *(cEV.getConnect(ic));
    char* eltTypConn = eltTypes[ic];
    nelts[ic] = cm.getSize();
    std::cout << "eltTypConn " << eltTypConn << std::endl;

    if (strcmp(eltTypConn, "BAR") == 0) 
    {
      nfaces[ic] = 2; nof[ic] = 1;
      face[ic].reserve(nfaces[ic] * nof[ic]);
      face[ic][0 + 0*nfaces[ic]] = 1;
      face[ic][1 + 0*nfaces[ic]] = 2;
    }
    else if (strcmp(eltTypConn, "TRI") == 0) 
    {
      nfaces[ic] = 3; nof[ic] = 2;
      face[ic].reserve(nfaces[ic] * nof[ic]);
      face[ic][0 + 0*nfaces[ic]] = 1; face[ic][0 + 1*nfaces[ic]] = 2;
      face[ic][1 + 0*nfaces[ic]] = 2; face[ic][1 + 1*nfaces[ic]] = 3;
      face[ic][2 + 0*nfaces[ic]] = 3; face[ic][2 + 1*nfaces[ic]] = 1;
    }
    else if (strcmp(eltTypConn, "QUAD") == 0) 
    {
      nfaces[ic] = 4; nof[ic] = 2;
      face[ic].reserve(nfaces[ic] * nof[ic]);
      face[ic][0 + 0*nfaces[ic]] = 1; face[ic][0 + 1*nfaces[ic]] = 2;
      face[ic][1 + 0*nfaces[ic]] = 2; face[ic][1 + 1*nfaces[ic]] = 3;
      face[ic][2 + 0*nfaces[ic]] = 3; face[ic][2 + 1*nfaces[ic]] = 4;
      face[ic][3 + 0*nfaces[ic]] = 4; face[ic][3 + 1*nfaces[ic]] = 1;
    }
    else if (strcmp(eltTypConn, "TETRA") == 0) 
    {
      nfaces[ic] = 4; nof[ic] = 3;
      face[ic].reserve(nfaces[ic] * nof[ic]);
      face[ic][0 + 0*nfaces[ic]] = 1; face[ic][0 + 1*nfaces[ic]] = 2; face[ic][0 + 2*nfaces[ic]] = 3;
      face[ic][1 + 0*nfaces[ic]] = 1; face[ic][1 + 1*nfaces[ic]] = 2; face[ic][1 + 2*nfaces[ic]] = 4;
      face[ic][2 + 0*nfaces[ic]] = 2; face[ic][2 + 1*nfaces[ic]] = 3; face[ic][2 + 2*nfaces[ic]] = 4;
      face[ic][3 + 0*nfaces[ic]] = 3; face[ic][3 + 1*nfaces[ic]] = 1; face[ic][3 + 2*nfaces[ic]] = 4;
    }
    else if (strcmp(eltTypConn, "PYRA") == 0) 
    {
      nfaces[ic] = 6; nof[ic] = 3; // 2 TRIs pour la base
      face[ic].reserve(nfaces[ic] * nof[ic]);
      face[ic][0 + 0*nfaces[ic]] = 1; face[ic][0 + 1*nfaces[ic]] = 4; face[ic][0 + 2*nfaces[ic]] = 3;
      face[ic][1 + 0*nfaces[ic]] = 3; face[ic][1 + 1*nfaces[ic]] = 2; face[ic][1 + 2*nfaces[ic]] = 1;
      face[ic][2 + 0*nfaces[ic]] = 1; face[ic][2 + 1*nfaces[ic]] = 2; face[ic][2 + 2*nfaces[ic]] = 5;
      face[ic][3 + 0*nfaces[ic]] = 2; face[ic][3 + 1*nfaces[ic]] = 3; face[ic][3 + 2*nfaces[ic]] = 5;
      face[ic][4 + 0*nfaces[ic]] = 3; face[ic][4 + 1*nfaces[ic]] = 4; face[ic][4 + 2*nfaces[ic]] = 5;
      face[ic][5 + 0*nfaces[ic]] = 4; face[ic][5 + 1*nfaces[ic]] = 1; face[ic][5 + 2*nfaces[ic]] = 5;
    }
    else if (strcmp(eltTypConn, "PENTA") == 0) 
    {
      nfaces[ic] = 5; nof[ic] = 4; // TRI degen
      face[ic].reserve(nfaces[ic] * nof[ic]);
      face[ic][0 + 0*nfaces[ic]] = 1; face[ic][0 + 1*nfaces[ic]] = 2; face[ic][0 + 2*nfaces[ic]] = 5; face[ic][0 + 3*nfaces[ic]] = 4;
      face[ic][1 + 0*nfaces[ic]] = 2; face[ic][1 + 1*nfaces[ic]] = 3; face[ic][1 + 2*nfaces[ic]] = 6; face[ic][1 + 3*nfaces[ic]] = 5;
      face[ic][2 + 0*nfaces[ic]] = 3; face[ic][2 + 1*nfaces[ic]] = 1; face[ic][2 + 2*nfaces[ic]] = 4; face[ic][2 + 3*nfaces[ic]] = 6;
      face[ic][3 + 0*nfaces[ic]] = 1; face[ic][3 + 1*nfaces[ic]] = 3; face[ic][3 + 2*nfaces[ic]] = 2; face[ic][3 + 3*nfaces[ic]] = 1;
      face[ic][4 + 0*nfaces[ic]] = 4; face[ic][4 + 1*nfaces[ic]] = 5; face[ic][4 + 2*nfaces[ic]] = 6; face[ic][4 + 3*nfaces[ic]] = 4;
    }
    else if (strcmp(eltTypConn, "HEXA") == 0) 
    {
      nfaces[ic] = 6; nof[ic] = 4;
      face[ic].reserve(nfaces[ic] * nof[ic]);
      face[ic][0 + 0*nfaces[ic]] = 1; face[ic][0 + 1*nfaces[ic]] = 4; face[ic][0 + 2*nfaces[ic]] = 3; face[ic][0 + 3*nfaces[ic]] = 2;
      face[ic][1 + 0*nfaces[ic]] = 1; face[ic][1 + 1*nfaces[ic]] = 2; face[ic][1 + 2*nfaces[ic]] = 6; face[ic][1 + 3*nfaces[ic]] = 5;
      face[ic][2 + 0*nfaces[ic]] = 2; face[ic][2 + 1*nfaces[ic]] = 3; face[ic][2 + 2*nfaces[ic]] = 7; face[ic][2 + 3*nfaces[ic]] = 6;
      face[ic][3 + 0*nfaces[ic]] = 3; face[ic][3 + 1*nfaces[ic]] = 4; face[ic][3 + 2*nfaces[ic]] = 8; face[ic][3 + 3*nfaces[ic]] = 7;
      face[ic][4 + 0*nfaces[ic]] = 1; face[ic][4 + 1*nfaces[ic]] = 5; face[ic][4 + 2*nfaces[ic]] = 8; face[ic][4 + 3*nfaces[ic]] = 4;
      face[ic][5 + 0*nfaces[ic]] = 5; face[ic][5 + 1*nfaces[ic]] = 6; face[ic][5 + 2*nfaces[ic]] = 7; face[ic][5 + 3*nfaces[ic]] = 8;
    }

    // Update total face count
    ntotfaces[ic+1] = ntotfaces[ic] + nelts[ic]*nfaces[ic];
  }

  for (size_t ic = 0; ic < eltTypes.size(); ic++) delete [] eltTypes[ic];

  E_Int nofmax = *max_element(begin(nof), end(nof));
  cFV.malloc(ntotfaces[nc], nofmax);

  // Boucle sur toutes les connectivites pour remplir cFV
  for (E_Int ic = 0; ic < nc; ic++)
  {
    FldArrayI& cm = *(cEV.getConnect(ic));
#pragma omp parallel for default (shared)
    for (E_Int e = 0; e < nelts[ic]; e++)
    {
      E_Int ind;
      for (E_Int f = 0; f < nfaces[ic]; f++)
      {
        ind = f + e*nfaces[ic] + ntotfaces[ic];
        for (E_Int i = 0; i < nof[ic]; i++) cFV(ind,i+1) = cm(e,face[ic][f + i*nfaces[ic]]);
      }
    }
  }
}