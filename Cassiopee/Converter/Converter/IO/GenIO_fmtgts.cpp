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

// Formated gts (Gnu Triangulated Surface lib) file support

# include <stdio.h>
# include <stdlib.h>
# include <vector>
# include <string.h>

# include "GenIO.h"
# include "Array/Array.h"
# include "Def/DefFunction.h"
# include "Connect/connect.h"
# include "CompGeom/compGeom.h"

using namespace K_FLD;
using namespace std;

//=============================================================================
/* gtsread 
   Read gts surfaces as BAR (edges) and TRI (surface) arrays. */
//=============================================================================
E_Int K_IO::GenIO::gtsread(
  char* file, char*& varString,
  vector<FldArrayF*>& structField,
  vector<E_Int>& ni, vector<E_Int>& nj, vector<E_Int>& nk,
  vector<FldArrayF*>& unstructField,
  vector<FldArrayI*>& connect,
  vector<E_Int>& eltType,
  vector<char*>& zoneNames)
{
  /* File Opening */
  FILE* ptrFile;
  ptrFile = fopen(file, "r");

  if (ptrFile == NULL)
  {
    printf("Warning: gtsread: cannot open file %s.\n", file);
    return 1;
  }

  varString = new char [8];
  strcpy(varString, "x,y,z");

  E_Int v = 1;

  while (v == 1) v = skipComment(ptrFile, '#');
  E_Int nv;
  readInt(ptrFile, nv);
  while (v == 1) v = skipComment(ptrFile, '#');
  E_Int nedges;
  readInt(ptrFile, nedges);
  while (v == 1) v = skipComment(ptrFile, '#');
  E_Int ne;
  readInt(ptrFile, ne, -1);
  skipLine(ptrFile);

  // vertices
  v = 1;
  while (v == 1) v = skipComment(ptrFile, '#');
  FldArrayF* coord = new FldArrayF(nv, 3);
  E_Float* x; E_Float* y; E_Float* z;
  x = coord->begin(1);
  y = coord->begin(2);
  z = coord->begin(3);
  for (E_Int i = 0; i < nv; i++)
  {
    readDouble(ptrFile, x[i]);
    readDouble(ptrFile, y[i]);
    readDouble(ptrFile, z[i]);
  }
  
  // edges
  v = 1;
  while (v == 1)
    v = skipComment(ptrFile, '#');
  FldArrayI edges(nedges, 2);
  E_Int* e1 = edges.begin(1);
  E_Int* e2 = edges.begin(2);
  for (E_Int i = 0; i < nedges; i++)
  {
    readInt(ptrFile, e1[i]);
    readInt(ptrFile, e2[i]);
  }

  // triangles
  while (v == 1) v = skipComment(ptrFile, '#');
  FldArrayI* c = new FldArrayI(ne, 3);
  E_Int* c1; E_Int* c2; E_Int* c3;
  c1 = c->begin(1);
  c2 = c->begin(2);
  c3 = c->begin(3);
  E_Int ind1, ind2, ind3;
  for (E_Int i = 0; i < ne; i++)
  {
    readInt(ptrFile, ind1);
    readInt(ptrFile, ind2);
    readInt(ptrFile, ind3);
    c1[i] = edges(ind1-1,1);
    c2[i] = edges(ind1-1,2);
    if (edges(ind2-1, 1) == c2[i])
      c3[i] = edges(ind2-1, 2);
    else if (edges(ind2-1, 1) == c1[i])
      c3[i] = edges(ind2-1, 2);
    else if (edges(ind2-1, 2) == c2[i])
      c3[i] = edges(ind2-1, 1);
    else 
      c3[i] = edges(ind2-1, 1);
  }
  unstructField.push_back(coord);
  connect.push_back(c);
  eltType.push_back(2);
  // nom de zone
  char* zoneName = new char [128];
  sprintf(zoneName, "Zone0");
  zoneNames.push_back(zoneName);
  fclose(ptrFile);
  return 0;
}

//=============================================================================
/* Write arrays as gts surface.
   Others are discarded. */
//=============================================================================
E_Int K_IO::GenIO::gtswrite(
  char* file, char* dataFmt, char* varString,
  vector<E_Int>& ni, vector<E_Int>& nj, vector<E_Int>& nk,
  vector<FldArrayF*>& structField,
  vector<FldArrayF*>& unstructField,
  vector<FldArrayI*>& connect,
  vector< vector<E_Int> >& eltTypes,
  vector<char*>& zoneNames)
{
  // All zones must have posx, posy, posz
  E_Int posx, posy, posz;
  posx = K_ARRAY::isCoordinateXPresent(varString);
  posy = K_ARRAY::isCoordinateYPresent(varString);
  posz = K_ARRAY::isCoordinateZPresent(varString);
  if (posx == -1 || posy == -1 || posz == -1)
  {
    printf("Warning: gtswrite: zones do not have coordinates. Not written.\n");
    return 1;
  }
  posx++; posy++; posz++;

  // Ecriture de l'entete
  FILE* ptrFile = fopen(file, "w");
  if (ptrFile == NULL) 
  {
    printf("gtswrite: I can't open file %s.\n", file);
    return 1;
  }

  // Calcul du tableau edge
  
  //fwrite(ptrFile, "%d %d %d\n", );
  printf("Warning: gtswrite: not implemented.\n");

  fclose(ptrFile);
  return 0;
}
