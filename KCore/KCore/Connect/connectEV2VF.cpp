/*    
    Copyright 2013-2023 Onera.

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
   Vertex->Faces. L'indice des faces est global, soit : nof + nelt*nfaces
   ou nfaces est le nbre de faces de l'elements, nelt le no de l'element
   (commencant a 0) et nof la numero local de la face 1,2,3,...
   cVF doit etre alloue a npts.
*/
//=============================================================================
void K_CONNECT::connectEV2VF(FldArrayI& cEV, const char* eltType,
                             vector< vector<E_Int> >& cVF)
{
  // Acces universel sur BE/ME
  E_Int nc = cEV.getNConnect();
  // Acces universel aux eltTypes
  vector<char*> eltTypes;
  K_ARRAY::extractVars(eltType, eltTypes);
  E_Int el_offset = 0; // element offset for subsequent connectivities
  // Number of elements and faces per connectivity
  vector<E_Int> nelts(nc);
  vector<E_Int> nfaces(nc);
  vector<E_Int> nof(nc);
  vector<vector<vector<E_Int> > > face;
  E_Int npts = cVF.size();

  // Size
  FldArrayI size(npts); size.setAllValuesAtNull();
  E_Int* sizep = size.begin();
  E_Int ind, indN;
    
  // Boucle sur toutes les connectivites pour remplir face et pre-evaluer
  // le nombre de faces connectees a chaque noeud
  for (E_Int ic = 0; ic < nc; ic++)
  {
    FldArrayI& cm = *(cEV.getConnect(ic));
    char* eltTypConn = eltTypes[ic];
    nelts[ic] = cm.getSize();

    vector<vector<E_Int> > facec(6);
    for (size_t ifc = 0; ifc < facec.size(); ifc++) facec[ifc].reserve(4);

    if (strcmp(eltTypConn, "BAR") == 0) 
    { 
      nfaces[ic] = 2; nof[ic] = 1;
      facec[0][0] = 1; facec[1][0] = 2;
    }
    else if (strcmp(eltTypConn, "QUAD") == 0) 
    {
      nfaces[ic] = 4; nof [ic]= 2;
      facec[0][0] = 1; facec[0][1] = 2;
      facec[1][0] = 2; facec[1][1] = 3;
      facec[2][0] = 3; facec[2][1] = 4;
      facec[3][0] = 4; facec[3][1] = 1;
    }
    else if (strcmp(eltTypConn, "TRI") == 0) 
    {
      nfaces[ic] = 3; nof[ic] = 2;
      facec[0][0] = 1; facec[0][1] = 2;
      facec[1][0] = 2; facec[1][1] = 3;
      facec[2][0] = 3; facec[2][1] = 1;
    }
    else if (strcmp(eltTypConn, "HEXA") == 0) 
    {
      nfaces[ic] = 6; nof[ic] = 4;
      facec[0][0] = 1; facec[0][1] = 4; facec[0][2] = 3; facec[0][3] = 2;
      facec[1][0] = 1; facec[1][1] = 2; facec[1][2] = 6; facec[1][3] = 5;
      facec[2][0] = 2; facec[2][1] = 3; facec[2][2] = 7; facec[2][3] = 6;
      facec[3][0] = 3; facec[3][1] = 4; facec[3][2] = 8; facec[3][3] = 7;
      facec[4][0] = 1; facec[4][1] = 5; facec[4][2] = 8; facec[4][3] = 4;
      facec[5][0] = 5; facec[5][1] = 6; facec[5][2] = 7; facec[5][3] = 8;
    }
    else if (strcmp(eltTypConn, "TETRA") == 0) 
    {
      nfaces[ic] = 4; nof[ic] = 3;
      facec[0][0] = 1; facec[0][1] = 3; facec[0][2] = 2;
      facec[1][0] = 1; facec[1][1] = 2; facec[1][2] = 4;
      facec[2][0] = 2; facec[2][1] = 3; facec[2][2] = 4;
      facec[3][0] = 3; facec[3][1] = 1; facec[3][2] = 4;
    }
    else if (strcmp(eltTypConn, "PYRA") == 0) 
    {
      nfaces[ic] = 5; nof[ic] = 3; // 2 TRIs pour la base
      facec[0][0] = 1; facec[0][1] = 4; facec[0][2] = 3;
      facec[1][0] = 3; facec[1][1] = 2; facec[1][2] = 1;
      facec[2][0] = 1; facec[2][1] = 2; facec[2][2] = 5; 
      facec[3][0] = 2; facec[3][1] = 3; facec[3][2] = 5;
      facec[4][0] = 3; facec[4][1] = 4; facec[4][2] = 5;
      facec[5][0] = 4; facec[5][1] = 1; facec[5][2] = 5;
    }
    else if (strcmp(eltTypConn, "PENTA") == 0) 
    {
      nfaces[ic] = 5; nof[ic] = 4; // TRI degen
      facec[0][0] = 1; facec[0][1] = 2; facec[0][2] = 5; facec[0][3] = 4;
      facec[1][0] = 2; facec[1][1] = 3; facec[1][2] = 6; facec[1][3] = 5;
      facec[2][0] = 3; facec[2][1] = 1; facec[2][2] = 4; facec[2][3] = 6;
      facec[3][0] = 1; facec[3][1] = 3; facec[3][2] = 2; facec[3][3] = 2;
      facec[4][0] = 4; facec[4][1] = 5; facec[4][2] = 6; facec[4][3] = 6;
    }

    face.push_back(facec);
  
    // Determine the number of faces connected to each vertex
    for (E_Int e = 0; e < nelts[ic]; e++)
    {
      for (E_Int f = 0; f < nfaces[ic]; f++)
      {
        for (E_Int i = 0; i < nof[ic]; i++) 
        {
          indN = cm(e,facec[f][i]);
          sizep[indN-1]++;
        }
      }
    }
  }

  // Boucle sur toutes les connectivites pour remplir cVF
  for (E_Int ic = 0; ic < nc; ic++)
  {
    FldArrayI& cm = *(cEV.getConnect(ic));
#pragma omp parallel for default (shared)
    for (E_Int i = 0; i < npts; i++) cVF[i].reserve(sizep[i]);  

    // cVF
    for (E_Int e = 0; e < nelts[ic]; e++)
    {
      for (E_Int f = 0; f < nfaces[ic]; f++)
      {
        ind = f + 1 + (el_offset + e)*nfaces[ic];
        for (E_Int i = 0; i < nof[ic]; i++) 
        {
          indN = cm(e,face[ic][f][i]);
          cVF[indN-1].push_back(ind);
        }
      }
    }
    el_offset += nelts[ic]; // increment element offset
  }
}