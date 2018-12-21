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
/* Change un connectivite Elts-Vertex (basic elements) en une connectivite
   Faces->Vertex. L'indice des faces est global, soit : nof + nelt*nfaces
   ou nfaces est le nbre de faces de l'elements, nelt le no de l'element
   (commencant a 0) et nof la numero local de la face 0,1,2,...
*/
//=============================================================================
void K_CONNECT::connectEV2FV(FldArrayI& cEV, char* eltType,
                             FldArrayI& cFV)
{
  E_Int face[6][4];
  E_Int nfaces = 0; E_Int nof = 0;
  if (strcmp(eltType, "BAR") == 0)
  { 
    nfaces = 2; nof = 1;
    face[0][0] = 1; face[1][0] = 2;
  }
  else if (strcmp(eltType, "QUAD") == 0) 
  {
    nfaces = 4; nof = 2;
    face[0][0] = 1; face[0][1] = 2;
    face[1][0] = 2; face[1][1] = 3;
    face[2][0] = 3; face[2][1] = 4;
    face[3][0] = 4; face[3][1] = 1;
  }
  else if (strcmp(eltType, "TRI") == 0) 
  {
    nfaces = 3; nof = 2;
    face[0][0] = 1; face[0][1] = 2;
    face[1][0] = 2; face[1][1] = 3;
    face[2][0] = 3; face[2][1] = 1;
  }
  else if (strcmp(eltType, "HEXA") == 0) 
  {
    nfaces = 6; nof = 4;
    face[0][0] = 1; face[0][1] = 4; face[0][2] = 3; face[0][3] = 2;
    face[1][0] = 1; face[1][1] = 2; face[1][2] = 6; face[1][3] = 5;
    face[2][0] = 2; face[2][1] = 3; face[2][2] = 7; face[2][3] = 6;
    face[3][0] = 3; face[3][1] = 4; face[3][2] = 8; face[3][3] = 7;
    face[4][0] = 1; face[4][1] = 5; face[4][2] = 8; face[4][3] = 4;
    face[5][0] = 5; face[5][1] = 6; face[5][2] = 7; face[5][3] = 8;
  }
  else if (strcmp(eltType, "TETRA") == 0) 
  {
    nfaces = 4; nof = 3;
    face[0][0] = 1; face[0][1] = 3; face[0][2] = 2;
    face[1][0] = 1; face[1][1] = 2; face[1][2] = 4;
    face[2][0] = 2; face[2][1] = 3; face[2][2] = 4;
    face[3][0] = 3; face[3][1] = 1; face[3][2] = 4;
  }
  else if (strcmp(eltType, "PYRA") == 0) 
  {
    nfaces = 5; nof = 3; // 2 TRIs pour la base
    face[0][0] = 1; face[0][1] = 4; face[0][2] = 3;
    face[1][0] = 3; face[1][1] = 2; face[1][2] = 1;
    face[2][0] = 1; face[2][1] = 2; face[2][2] = 5; 
    face[3][0] = 2; face[3][1] = 3; face[3][2] = 5;
    face[4][0] = 3; face[4][1] = 4; face[4][2] = 5;
    face[5][0] = 4; face[5][1] = 1; face[5][2] = 5;
  }
  else if (strcmp(eltType, "PENTA") == 0) 
  {
    nfaces = 5; nof = 4; // TRI degen
    face[0][0] = 1; face[0][1] = 2; face[0][2] = 5; face[0][3] = 4;
    face[1][0] = 2; face[1][1] = 3; face[1][2] = 6; face[1][3] = 5;
    face[2][0] = 3; face[2][1] = 1; face[2][2] = 4; face[2][3] = 6;
    face[3][0] = 1; face[3][1] = 3; face[3][2] = 2; face[3][3] = 2;
    face[4][0] = 4; face[4][1] = 5; face[4][2] = 6; face[4][3] = 6;
  }

  E_Int nelts = cEV.getSize();
  E_Int nftot = nelts * nfaces;
  cFV.malloc(nftot, nof);
  
#pragma omp parallel for default (shared)
  for (E_Int e = 0; e < nelts; e++)
  {
    E_Int ind;
    for (E_Int f = 0; f < nfaces; f++)
    {
      ind = f + e*nfaces;
      for (E_Int i = 0; i < nof; i++) cFV(ind,i+1) = cEV(e,face[f][i]);
    }
  }
}
